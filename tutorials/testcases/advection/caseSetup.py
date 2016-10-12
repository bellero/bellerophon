#!/usr/bin/env python
  
import shutil
import os

import helperFuncs

# user input:

class Setup:
  def __init__(self):
    # Number of cells in vertical direction
    self.nys=[10, 20, 50, 75, 100, 150, 200]
    #self.nys=[200]
  
    # Interpolation methods
    #self.interpolationMethods=["donorCellTet"]
    self.interpolationMethods=["donorCell", "donorCellGradient", "donorCellTet"]
    #self.interpolationMethods=["donorCellGradient", "donorCellTet"]
  
    # Divergence schemes
    #self.divSchemes=["limitedCubic01 0.0 1.0",
    #                 "SFCD",
    #                 "filteredLinear",
    #                 "vanLeer",
    #                 "upwind",
    #                 "linear",
    #                 "limitedLinear 1",
    #                 "linearUpwind Gauss linear"]
    self.divSchemes=["limitedCubic01 0.0 1.0", "limitedLinear 1"]
    #self.divSchemes=["limitedLinear 1"]

    # Length of the domain
    self.length=2.5
  
    # Exact mass from bounding box
    self.massE=0.012
  
    self.rootDir=os.getcwd()+"/"
  
    # Name of template directory
    self.templateDir=self.rootDir+"templates/"
  
    # Name of calculations directory
    self.calculationsDir=self.rootDir+"calculations/"
  
    # Name of postprocessing directory
    self.postDir=self.rootDir+"post/"
  
    # Name for results
    self.resultDir=self.postDir+"results/"
  
    self.latestTime=10000
    self.nProcs=5

def doSingle(s):
  """
  @type s: Setup
  """
  for ny in s.nys:
    simpleCaseDir=s.calculationsDir+"single/simple/N"+str(ny)+"/"
    helperFuncs.copyTemplate(s.templateDir+"simple/single",simpleCaseDir)
    helperFuncs.createSingleBlockMesh(simpleCaseDir, ny, s.length)
  
    os.chdir(simpleCaseDir)
    helperFuncs.runApp("decomposePar")
    helperFuncs.runAppParallel("simpleFoam",s.nProcs,"-parallel")
    helperFuncs.runApp("reconstructPar")
    helperFuncs.runApp("sample","-latestTime")
    shutil.move("postProcessing/sets",s.resultDir+"U_single_N"+str(ny))
    
    shutil.rmtree(str(s.latestTime)+"/uniform")
  
    for scheme in s.divSchemes:

      transportCaseDir=s.calculationsDir+"single/transport/"+helperFuncs.caseName("",ny,scheme)+"/"
      
      helperFuncs.copyTemplate(s.templateDir+"transport/",transportCaseDir)
      os.chdir(transportCaseDir)
      os.system("cp "+simpleCaseDir+"constant/polyMesh/* "+transportCaseDir+"constant/polyMesh")
      os.system("cp "+simpleCaseDir+str(s.latestTime)+"/* "+transportCaseDir+"0")
      os.system("sed -i 's/DIVSCHEME/"+scheme+"/' system/fvSchemes")
      os.system("sed -i 's/GRADSCHEME/Gauss linear/' system/fvSchemes")
      helperFuncs.runApp("setFields")
      os.system("cp system/single/* system/")
      helperFuncs.runApp("decomposePar")
      helperFuncs.runAppParallel("scalarTransportFoam",s.nProcs,"-parallel")
      helperFuncs.runApp("reconstructPar")
      helperFuncs.runApp("sample","-latestTime")

      resultName="_".join((s.resultDir+"single_scalarDist_"+scheme+"_N"+str(ny)+".dat").split())
      shutil.move("postProcessing/sets/1.5/sampleCloud_T.xy",resultName)

      os.chdir(s.rootDir)
  print("Finished single")


def doPerfectMatching(s):
  """
  @type s Setup
  """
  for ny in s.nys:
    for method in s.interpolationMethods:
      simpleCaseDir=s.calculationsDir+"perfect/simple/"+helperFuncs.caseName(method,ny)+"/"
      helperFuncs.copyTemplate(s.templateDir+"simple/structured",simpleCaseDir)
      helperFuncs.createPerfectMatchBlockMesh(simpleCaseDir, ny, s.length)
    
      os.chdir(simpleCaseDir)
      os.system("sed -i 's/METHOD/"+method+"/' constant/bellerophonDict")
      helperFuncs.runApp("decomposePar")
      helperFuncs.runAppParallel("simpleFoam",s.nProcs,"-parallel")
      helperFuncs.runApp("reconstructPar")
      helperFuncs.runApp("sample","-latestTime")
      shutil.move("postProcessing/sets",s.resultDir+"U_perfect_"+helperFuncs.caseName(method,ny))

      resultName="_".join((s.resultDir+"perfect_flux_"+method+"_N"+str(ny)+".dat").split())
      helperFuncs.extractFluxConservation(simpleCaseDir,"log.simpleFoam",resultName)
      
      shutil.rmtree(str(s.latestTime)+"/uniform")
    
      for scheme in s.divSchemes:
  
        transportCaseDir=s.calculationsDir+"perfect/transport/"+helperFuncs.caseName(method,ny,scheme)+"/"
        
        helperFuncs.copyTemplate(s.templateDir+"transport/",transportCaseDir)
        os.chdir(transportCaseDir)
        os.system("cp "+simpleCaseDir+"constant/polyMesh/* "+transportCaseDir+"constant/polyMesh")
        os.system("cp "+simpleCaseDir+str(s.latestTime)+"/* "+transportCaseDir+"0")
        os.system("sed -i 's/METHOD/"+method+"/' constant/bellerophonDict")
        os.system("sed -i 's/DIVSCHEME/"+scheme+"/' system/fvSchemes")
        os.system("sed -i 's/GRADSCHEME/Gauss linear/' system/fvSchemes")
        helperFuncs.runApp("setFields")
        helperFuncs.runApp("decomposePar")
        helperFuncs.runAppParallel("scalarTransportFoam",s.nProcs,"-parallel")
        helperFuncs.runApp("reconstructPar")
        helperFuncs.runApp("sample","-latestTime")
  
        resultName="_".join((s.resultDir+"perfect_scalarDist_"+scheme+"_"+method+"_N"+str(ny)+".dat").split())
        shutil.move("postProcessing/sets/1.5/sampleCloud_T.xy",resultName)
        resultName="_".join((s.resultDir+"perfect_scalar_"+scheme+"_"+method+"_N"+str(ny)+".dat").split())
        helperFuncs.extractScalarConservation(transportCaseDir,resultName, s.massE,ny)
        
        os.chdir(s.rootDir)
  print("Finished perfect")


def doShifted(s):
  ## Next: shifted mesh
  for ny in s.nys:
    for method in s.interpolationMethods:
      simpleCaseDir=s.calculationsDir+"shifted/simple/"+helperFuncs.caseName(method,ny)+"/"
      helperFuncs.copyTemplate(s.templateDir+"simple/structured",simpleCaseDir)
      helperFuncs.createShiftedBlockMesh(simpleCaseDir, ny, s.length)
    
      os.chdir(simpleCaseDir)
      os.system("sed -i 's/METHOD/"+method+"/' constant/bellerophonDict")
      helperFuncs.runApp("decomposePar")
      helperFuncs.runAppParallel("simpleFoam",s.nProcs,"-parallel")
      helperFuncs.runApp("reconstructPar")
      helperFuncs.runApp("sample","-latestTime")
      shutil.move("postProcessing/sets",s.resultDir+"U_shifted_"+helperFuncs.caseName(method,ny))

      resultName="_".join((s.resultDir+"shifted_flux_"+"_"+method+"_N"+str(ny)+".dat").split())
      helperFuncs.extractFluxConservation(simpleCaseDir,"log.simpleFoam",resultName)
      
      shutil.rmtree(str(s.latestTime)+"/uniform")
    
      for scheme in s.divSchemes:
  
        transportCaseDir=s.calculationsDir+"shifted/transport/"+helperFuncs.caseName(method,ny,scheme)+"/"
        
        helperFuncs.copyTemplate(s.templateDir+"transport/",transportCaseDir)
        os.chdir(transportCaseDir)
        os.system("cp "+simpleCaseDir+"constant/polyMesh/* "+transportCaseDir+"constant/polyMesh")
        os.system("cp "+simpleCaseDir+str(s.latestTime)+"/* "+transportCaseDir+"0")
        os.system("sed -i 's/METHOD/"+method+"/' constant/bellerophonDict")
        os.system("sed -i 's/DIVSCHEME/"+scheme+"/' system/fvSchemes")
        os.system("sed -i 's/GRADSCHEME/Gauss linear/' system/fvSchemes")
        helperFuncs.runApp("setFields")
        helperFuncs.runApp("decomposePar")
        helperFuncs.runAppParallel("scalarTransportFoam",s.nProcs,"-parallel")
        helperFuncs.runApp("reconstructPar")
        helperFuncs.runApp("sample","-latestTime")
  
        resultName="_".join((s.resultDir+"shifted_scalarDist_"+scheme+"_"+method+"_N"+str(ny)+".dat").split())
        shutil.move("postProcessing/sets/1.5/sampleCloud_T.xy",resultName)
        resultName="_".join((s.resultDir+"shifted_scalar_"+scheme+"_"+method+"_N"+str(ny)+".dat").split())
        helperFuncs.extractScalarConservation(transportCaseDir,resultName, s.massE,ny)
        
        os.chdir(s.rootDir)
  print("Finished shifted")

def doTri(s):  
  ## unstructured tri
  for ny in s.nys:
    for method in s.interpolationMethods:
      simpleCaseDir=s.calculationsDir+"tri/simple/"+helperFuncs.caseName(method,ny)+"/"
      helperFuncs.copyTemplate(s.templateDir+"simple/unstructured",simpleCaseDir)
      helperFuncs.createUnstructured3DTriMeshes(simpleCaseDir, ny, s.length)
      os.chdir(simpleCaseDir)
  
      os.system("sed -i 's/METHOD/"+method+"/' constant/bellerophonDict")
      helperFuncs.runApp("extrudeMesh")
      helperFuncs.runApp("changeDictionary")
      helperFuncs.runApp("topoSet")
      helperFuncs.runApp("holeCutter")
      helperFuncs.runApp("subsetMesh","liveCells -overwrite")
      helperFuncs.runApp("createPatch","-overwrite")
      if os.path.isdir(simpleCaseDir+"0"):
        shutil.rmtree(simpleCaseDir+"0/")
      shutil.copytree(simpleCaseDir+"0.org",simpleCaseDir+"0")
      helperFuncs.runApp("decomposePar")
      helperFuncs.runAppParallel("simpleFoam",s.nProcs,"-parallel")
      helperFuncs.runApp("reconstructPar")
      helperFuncs.runApp("sample","-latestTime")
      shutil.move("postProcessing/sets",s.resultDir+"U_tri_"+helperFuncs.caseName(method,ny))

      resultName="_".join((s.resultDir+"tri_flux_"+"_"+method+"_N"+str(ny)+".dat").split())
      helperFuncs.extractFluxConservation(simpleCaseDir,"log.simpleFoam",resultName)
      
      shutil.rmtree(str(s.latestTime)+"/uniform")
      
      for scheme in s.divSchemes:
  
        transportCaseDir=s.calculationsDir+"tri/transport/"+helperFuncs.caseName(method,ny,scheme)+"/"
        
        helperFuncs.copyTemplate(s.templateDir+"transport/",transportCaseDir)
        os.chdir(transportCaseDir)
        os.system("cp "+simpleCaseDir+"constant/polyMesh/* "+transportCaseDir+"constant/polyMesh")
        os.system("cp "+simpleCaseDir+str(s.latestTime)+"/* "+transportCaseDir+"0")
        os.system("sed -i 's/METHOD/"+method+"/' constant/bellerophonDict")
        os.system("sed -i 's/DIVSCHEME/"+scheme+"/' system/fvSchemes")
        os.system("sed -i 's/GRADSCHEME/leastSquares/' system/fvSchemes")
        helperFuncs.runApp("setFields")
        helperFuncs.runApp("decomposePar")
        helperFuncs.runAppParallel("scalarTransportFoam",s.nProcs,"-parallel")
        helperFuncs.runApp("reconstructPar")
        helperFuncs.runApp("sample","-latestTime")

        resultName="_".join((s.resultDir+"tri_scalarDist_"+scheme+"_"+method+"_N"+str(ny)+".dat").split())
        shutil.move("postProcessing/sets/1.5/sampleCloud_T.xy",resultName)
        resultName="_".join((s.resultDir+"tri_scalar_"+scheme+"_"+method+"_N"+str(ny)+".dat").split())
        helperFuncs.extractScalarConservation(transportCaseDir,resultName, s.massE,ny)
        
        os.chdir(s.rootDir)
  print("Finished tri")

def doPoly(s):  
  ## Finally unstructured poly
  for ny in s.nys:
    for method in s.interpolationMethods:
      simpleCaseDir=s.calculationsDir+"poly/simple/"+helperFuncs.caseName(method,ny)+"/"
      helperFuncs.copyTemplate(s.templateDir+"simple/unstructured",simpleCaseDir)
      helperFuncs.createUnstructured3DTriMeshes(simpleCaseDir, ny, s.length)
      os.chdir(simpleCaseDir)
  
      os.system("sed -i 's/METHOD/"+method+"/' constant/bellerophonDict")
      helperFuncs.runApp("polyDualMesh","89 -overwrite -concaveMultiCells")
      helperFuncs.runApp("extrudeMesh")
      os.system("rm -r constant/polyMesh/sets constant/polyMesh/*Zones*")
      helperFuncs.runApp("splitMeshRegions","-makeCellZones -overwrite")
      os.system("rm -r constant/polyMesh/*Zones*")
      helperFuncs.runApp("topoSet","-dict system/topoSetDict.1")
      helperFuncs.runApp("changeDictionary")
      helperFuncs.runApp("topoSet")
      helperFuncs.runApp("holeCutter")
      if os.path.isdir(simpleCaseDir+"0"):
        shutil.rmtree(simpleCaseDir+"0/")
        os.mkdir(simpleCaseDir+"0")
      helperFuncs.runApp("subsetMesh","liveCells -overwrite")
      helperFuncs.runApp("createPatch","-overwrite")
      if os.path.isdir(simpleCaseDir+"0"):
        shutil.rmtree(simpleCaseDir+"0/")
      shutil.copytree(simpleCaseDir+"0.org",simpleCaseDir+"0")
      helperFuncs.runApp("decomposePar")
      helperFuncs.runAppParallel("simpleFoam",s.nProcs,"-parallel")
      helperFuncs.runApp("reconstructPar")
      helperFuncs.runApp("sample","-latestTime")
      shutil.move("postProcessing/sets",s.resultDir+"U_poly_"+helperFuncs.caseName(method,ny))
      
      resultName="_".join((s.resultDir+"poly_flux_"+"_"+method+"_N"+str(ny)+".dat").split())
      helperFuncs.extractFluxConservation(simpleCaseDir,"log.simpleFoam",resultName)
      
      shutil.rmtree(str(s.latestTime)+"/uniform")

      for scheme in s.divSchemes:
  
        transportCaseDir=s.calculationsDir+"poly/transport/"+helperFuncs.caseName(method,ny,scheme)+"/"
        
        helperFuncs.copyTemplate(s.templateDir+"transport/",transportCaseDir)
        os.chdir(transportCaseDir)
        os.system("cp "+simpleCaseDir+"constant/polyMesh/* "+transportCaseDir+"constant/polyMesh")
        os.system("cp "+simpleCaseDir+str(s.latestTime)+"/* "+transportCaseDir+"0")
        os.system("sed -i 's/METHOD/"+method+"/' constant/bellerophonDict")
        os.system("sed -i 's/DIVSCHEME/"+scheme+"/' system/fvSchemes")
        os.system("sed -i 's/GRADSCHEME/leastSquares/' system/fvSchemes")
        helperFuncs.runApp("setFields")
        helperFuncs.runApp("decomposePar")
        helperFuncs.runAppParallel("scalarTransportFoam",s.nProcs,"-parallel")
        helperFuncs.runApp("reconstructPar")
        helperFuncs.runApp("sample","-latestTime")

        resultName="_".join((s.resultDir+"poly_scalarDist_"+scheme+"_"+method+"_N"+str(ny)+".dat").split())
        shutil.move("postProcessing/sets/1.5/sampleCloud_T.xy",resultName)
        resultName="_".join((s.resultDir+"poly_scalar_"+scheme+"_"+method+"_N"+str(ny)+".dat").split())
        helperFuncs.extractScalarConservation(transportCaseDir,resultName, s.massE,ny)
        
        os.chdir(s.rootDir)
  print("Finished poly")
