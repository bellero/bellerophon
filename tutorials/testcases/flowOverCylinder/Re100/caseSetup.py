#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue May  5 13:38:49 2015

@author: agross

Two-dimensional, laminar flow around round cylinder using overset and
consistent grid.

Performs a validation test by grid refinement and compares values for drag and
lift in average as well as Strouhal number with experimental results
"""

import math
import os
import shutil
import prepostprocessing


class Setup:
    def __init__(self):
        self.nus = [20, 40, 60, 70] # Cells per forth of the diameter
        #self.nus = [20, 40]

        self.methods = ["donorCell", "donorCellGradient", "donorCellTet"]
        #self.methods = ["donorCellTet"]

        self.r0 = 0.5  # Radius of cylinder
        self.r1 = 2.5  # Radius of cylinder domain
        self.r2 = 20.0  # Radius of overall domain
        self.AR = 2.0  # Aspect ratio of cells on the surface of the cylinder

        self.rootPath = os.getcwd() + "/"
        self.calculations = self.rootPath + "calculations/"
        self.templates = self.rootPath + "templates/"
        self.results = self.rootPath + "post/results/"
        self.nProcs = 3


def doSingle(s):
    print("Running single cases")

    for idx, nu in enumerate(s.nus):

        print("  Running single case N" + str(nu))
        deltaU0 = 0.5 * math.pi * s.r0 / (s.AR * nu)  # Height of cells on cylinder surface
        deltaU1 = 0.5 * math.pi * s.r1 / nu  # Height of cells on block boundary
        deltaU2 = 0.5 * math.pi * s.r2 / nu  # Height of cells on outer cylinder

        print("deltaU0 : " + str(deltaU0))
        print("deltaU1 : " + str(deltaU1))
        print("deltaU2 : " + str(deltaU2))

        E1 = deltaU1 / deltaU0
        E2 = deltaU2 / deltaU1

        print("E1      : " + str(E1))
        print("E2      : " + str(E2))

        nr1 = prepostprocessing.nFromLEDelta(s.r1 - s.r0, E1, deltaU0)  # Number of cells in radial direction (inner)
        nr2 = prepostprocessing.nFromLEDelta(s.r2 - s.r1, E2, deltaU1)  # Number of cells in radial direction (outer)

        print("nr1     : " + str(nr1))
        print("nr2     : " + str(nr2))

        transientCaseDir = s.calculations + "single_N" + str(nu) + "_transient"
        shutil.copytree(s.templates + "single_transient", transientCaseDir)

        if idx is 0:
            steadyCaseDir = s.calculations + "single_N" + str(nu) + "_steady"
            shutil.copytree(s.templates + "single_steady", steadyCaseDir)
            prepostprocessing.createSingleBlockMesh(steadyCaseDir, nu, nr1, nr2, E1, E2, s.r0, s.r1, s.r2)
            os.chdir(steadyCaseDir)
            prepostprocessing.decompAndRun("simpleFoam", s.nProcs)
            if os.path.isdir("500/uniform"):
                shutil.rmtree("500/uniform")

            os.chdir(transientCaseDir)
            shutil.copytree(steadyCaseDir + "/500", transientCaseDir + "/0")
            shutil.copytree(steadyCaseDir + "/constant/polyMesh", transientCaseDir + "/constant/polyMesh")
        else:
            prepostprocessing.createSingleBlockMesh(transientCaseDir, nu, nr1, nr2, E1, E2, s.r0, s.r1, s.r2)
            os.chdir(transientCaseDir)
            shutil.copytree(transientCaseDir + "/0.org", transientCaseDir + "/0")
            prepostprocessing.runApp("mapFields",
                                     "-sourceTime 500 -fields '(U p)' -mapMethod mapNearest " + s.calculations + "single_N" + str(
                                         s.nus[idx - 1]) + "_transient")

        prepostprocessing.decompAndRun("pimpleFoam", s.nProcs)
        shutil.copyfile("postProcessing/forces/0/forceCoeffs.dat", s.results + "forces_single_N" + str(nu) + ".dat")
        print("  Finished single N" + str(nu))
    print("Finished single")


def doOverset(s):
    print("Running overset cases")
    for idx, nu in enumerate(s.nus):

        print("  Running overset case N" + str(nu))
        deltaU0 = 0.5 * math.pi * s.r0 / (s.AR * nu)  # Height of cells on cylinder surface
        deltaU1 = 0.5 * math.pi * s.r1 / nu  # Height of cells on block boundary

        E = deltaU1 / deltaU0
        nr = prepostprocessing.nFromLEDelta(s.r1 - s.r0, E, deltaU0)  # Number of cells in radial direction (inner)
        nb = int(0.8 * s.r2 / deltaU1)  # Number of cells aloung background domain

        for method in s.methods:
            print("  Running overset case N" + str(nu) + "_" + method)
            transientCaseDir = s.calculations + "overset_N" + str(nu) + "_" + method + "_transient"
            shutil.copytree(s.templates + "overset_transient", transientCaseDir)

            if idx is 0:
                steadyCaseDir = s.calculations + "overset_N" + str(nu) + "_" + method + "_steady"
                shutil.copytree(s.templates + "overset_steady", steadyCaseDir)
                prepostprocessing.createOversetBlockMesh(steadyCaseDir, nu, nr, nb, E, s.r0, s.r1, s.r2)
                os.chdir(steadyCaseDir)
                os.mkdir("0")
                os.system("sed -i \"s/METHOD/" + method + "/\" constant/bellerophonDict")
                prepostprocessing.runApp("setSet", "-batch setBatch0")
                prepostprocessing.runApp("refineMesh", "-dict system/refineMeshDict -overwrite")
                os.system("rm -rf constant/polyMesh/sets")
                prepostprocessing.runApp("setSet", "-batch setBatch1")
                prepostprocessing.runApp("refineMesh", "-dict system/refineMeshDict -overwrite")
                os.system("rm -rf constant/polyMesh/sets")
                prepostprocessing.runApp("setSet", "-batch setBatch2")
                os.system("rm -rf constant/polyMesh/sets")
                prepostprocessing.runApp("surfaceMeshTriangulate", "-patches \"(inner)\" inner.stl")
                prepostprocessing.runApp("holeCutter", "-cellState")
                prepostprocessing.runApp("subsetMesh", "-overwrite liveCells")
                prepostprocessing.runApp("createPatch", "-overwrite")
                if os.path.isdir("0"):
                    shutil.rmtree("0")
                shutil.copytree("0.org", "0")
                prepostprocessing.decompAndRun("simpleFoam", s.nProcs)
                if os.path.isdir("500/uniform"):
                    shutil.rmtree("500/uniform")

                os.chdir(transientCaseDir)
                shutil.copytree(steadyCaseDir + "/500", transientCaseDir + "/0")
                shutil.copytree(steadyCaseDir + "/constant/polyMesh", transientCaseDir + "/constant/polyMesh")
            else:
                prepostprocessing.createOversetBlockMesh(transientCaseDir, nu, nr, nb, E, s.r0, s.r1, s.r2)
                os.chdir(transientCaseDir)
                os.mkdir("0")
                os.system("sed -i \"s/METHOD/" + method + "/\" constant/bellerophonDict")
                prepostprocessing.runApp("setSet", "-batch setBatch0")
                prepostprocessing.runApp("refineMesh", "-dict system/refineMeshDict -overwrite")
                os.system("rm -rf constant/polyMesh/sets")
                prepostprocessing.runApp("setSet", "-batch setBatch1")
                prepostprocessing.runApp("refineMesh", "-dict system/refineMeshDict -overwrite")
                os.system("rm -rf constant/polyMesh/sets")
                prepostprocessing.runApp("setSet", "-batch setBatch2")
                os.system("rm -rf constant/polyMesh/sets")
                prepostprocessing.runApp("surfaceMeshTriangulate", "-patches \"(inner)\" inner.stl")
                prepostprocessing.runApp("holeCutter", "-cellState")
                prepostprocessing.runApp("subsetMesh", "-overwrite liveCells")
                prepostprocessing.runApp("createPatch", "-overwrite")
                if os.path.isdir("0"):
                    shutil.rmtree("0")
                shutil.copytree("0.org", "0")
                prepostprocessing.runApp("mapFields",
                                         "-sourceTime 500 -fields '(U p)' -mapMethod mapNearest " + s.calculations + "overset_N" + str(
                                             s.nus[idx - 1]) + "_" + method + "_transient")

            os.system("sed -i \"s/METHOD/" + method + "/\" constant/bellerophonDict")
            prepostprocessing.decompAndRun("pimpleFoam", s.nProcs)
            prepostprocessing.extractFluxConservation(transientCaseDir, "log.pimpleFoam",
                                                      s.results + "conservation_overset_" + method + "_N" + str(
                                                          nu) + ".dat")
            shutil.copyfile("postProcessing/forces/0/forceCoeffs.dat",
                            s.results + "forces_overset_" + method + "_N" + str(nu) + ".dat")
        print("  Finished overset N" + str(nu))
    print("Finished overset")


def doPostprocessing(s):
    for nu in s.nus:
        prepostprocessing.postprocessForces(s.results + "forces_single_N" + str(nu) + ".dat", s.results, "single", nu)
        for method in s.methods:
            prepostprocessing.postprocessForces(s.results + "forces_overset_" + method + "_N" + str(nu) + ".dat",
                                                s.results, method, nu)
      
