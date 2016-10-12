# -*- coding: utf-8 -*-
"""
Created on Wed May  6 22:27:21 2015

@author: andreas
"""

import numpy as np
import math
import os

def postprocessForces(forceFile,postPath,method,n):
  """
  read Force coeffs file and extract data:
    * pdf
    * strouhal number
    * average CD
    * last period
  """
  inFile = open(forceFile,'r')
  t=[]
  CD=[]
  CL=[]
  for l in inFile:
    if l[0] != "#":
      line = map(float,l.split())
      t.append(line[0])
      CD.append(line[2])
      CL.append(line[3])
  inFile.close()

  nData = int(len(t)/2)
  t  = t[nData:]
  CD = CD[nData:]
  CL = CL[nData:]

  avCD = np.average(CD)
  avCL = np.average(CL)

  upsampledT =np.linspace(t[0],t[-1],len(t)*10)
  upsampledCD=np.interp(upsampledT,t,[c-avCD for c in CD])
 
  Ns = len(upsampledT)
  DeltaT = upsampledT[1]-upsampledT[0]
  spectrumCD = 2.0/Ns*np.abs(np.fft.fft(upsampledCD))
  freqs = np.fft.fftfreq(Ns,DeltaT)
  idx = np.argsort(freqs)

  spectrumCD = spectrumCD[idx][int(Ns/2):]
  freqs = freqs[idx][int(Ns/2):]
  
  Sr = 0.5*freqs[np.argmax(spectrumCD)]
  
  print("Strouhal number is "+str(Sr))
  
  ## Output
  path = postPath+"/"+method+"_N"+str(n)+"_"
  spectraFile = open(path+"spectra.dat",'w')
  spectraFile.write("#f[hz]\tCD\n")
  for idx in xrange(len(freqs)):
    spectraFile.write(str(freqs[idx])+"\t"+str(spectrumCD[idx])+"\n")
  spectraFile.close()
  
  lastPeriodFile = open(path+"lastPeriod.dat",'w')
  T=1.0/Sr
  T0=t[-1]-T
  lastPeriodFile.write("#t\tCD\tCL\n")
  for idx, ti in enumerate(t):
    if ti > T0:
      lastPeriodFile.write(str((ti-T0)/T)+"\t"+str(CD[idx])+"\t"+str(CL[idx])+"\n")
  lastPeriodFile.close()
  
  path = postPath+"/"+method+"_"
  SrFile = open(path+"Sr.dat",'a')
  SrFile.write(str(n)+"\t"+str(Sr)+"\n")
  SrFile.close()  
  
  CDFile = open(path+"CD.dat",'a')
  CDFile.write(str(n)+"\t"+str(avCD)+"\n")

  CLFile = open(path+"CL.dat",'a')
  CLFile.write(str(n)+"\t"+str(avCL)+"\n")

def LFromNDeltaEpsilon(N, delta, epsilon):
  """
  Calculate edge length from number of cells, cell size and expansion ratio
  """
  return delta*(1.0-math.pow(epsilon,N))/(1.0-epsilon)

def nFromLEpsilonDelta(L, epsilon, Delta):
  """
  Calculate number of cells along edge from length, expansion ratio and size
  of first cell
  """
  return int(math.ceil((math.log((epsilon-1.0)*L/Delta+1.0)/math.log(epsilon))))
  
def nFromLEDelta(L, E, Delta):
  """
  Calculate number of cells along edge from length, ratio between first and
  last cell and Size of first cell
  """
  epsilon = 0.0
  n = int(math.ceil(L/Delta))
  while n > 2 and epsilon < 1.0:
    eps = math.pow(E,1.0/float(n-1))
    length2=LFromNDeltaEpsilon(n,Delta,eps)
    if(L > length2):
      n=n+1
      epsilon=math.pow(E,1.0/float(n-1))
    else:
      n=n-1
  return int(n)

def createSingleBlockMesh(path, nu, nr1, nr2, E1, E2, r0, r1, r2):
  """
  Create blockMesh with two concentric regions
  """
  if not os.path.exists(path+"/constant/polyMesh"):
    os.makedirs(path+"/constant/polyMesh")
  
  c0=math.sqrt(0.5)*r0
  c1=math.sqrt(0.5)*r1
  c2=math.sqrt(0.5)*r2
  height=1.0

  blockMeshFile = open(path+"/constant/polyMesh/blockMeshDict",'w')
  blockMeshFile.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
  blockMeshFile.write("| =========                 |                                                 |\n")
  blockMeshFile.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
  blockMeshFile.write("|  \\    /   O peration     | Version:  2.1.x                                 |\n")
  blockMeshFile.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")
  blockMeshFile.write("|    \\/     M anipulation  |                                                 |\n")
  blockMeshFile.write("\*---------------------------------------------------------------------------*/\n")
  blockMeshFile.write("FoamFile\n")
  blockMeshFile.write("{\n")
  blockMeshFile.write("    version     2.0;\n")
  blockMeshFile.write("    format      ascii;\n")
  blockMeshFile.write("    class       dictionary;\n")
  blockMeshFile.write("    object      blockMeshDict;\n")
  blockMeshFile.write("}\n")
  blockMeshFile.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("convertToMeters 1;\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("vertices\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write("    ("+str( c0)+" "+str( c0)+" "+str( 0.0  )+") // 0\n")
  blockMeshFile.write("    ("+str(-c0)+" "+str( c0)+" "+str( 0.0  )+") // 1\n")
  blockMeshFile.write("    ("+str(-c0)+" "+str(-c0)+" "+str( 0.0  )+") // 2\n")
  blockMeshFile.write("    ("+str( c0)+" "+str(-c0)+" "+str( 0.0  )+") // 3\n")
  blockMeshFile.write("    ("+str( c1)+" "+str( c1)+" "+str( 0.0  )+") // 4\n")
  blockMeshFile.write("    ("+str(-c1)+" "+str( c1)+" "+str( 0.0  )+") // 5\n")
  blockMeshFile.write("    ("+str(-c1)+" "+str(-c1)+" "+str( 0.0  )+") // 6\n")
  blockMeshFile.write("    ("+str( c1)+" "+str(-c1)+" "+str( 0.0  )+") // 7\n")
  blockMeshFile.write("    ("+str( c2)+" "+str( c2)+" "+str( 0.0  )+") // 8\n")
  blockMeshFile.write("    ("+str(-c2)+" "+str( c2)+" "+str( 0.0  )+") // 9\n")
  blockMeshFile.write("    ("+str(-c2)+" "+str(-c2)+" "+str( 0.0  )+") //10\n")
  blockMeshFile.write("    ("+str( c2)+" "+str(-c2)+" "+str( 0.0  )+") //11\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    ("+str( c0)+" "+str( c0)+" "+str(height)+") //12\n")
  blockMeshFile.write("    ("+str(-c0)+" "+str( c0)+" "+str(height)+") //13\n")
  blockMeshFile.write("    ("+str(-c0)+" "+str(-c0)+" "+str(height)+") //14\n")
  blockMeshFile.write("    ("+str( c0)+" "+str(-c0)+" "+str(height)+") //15\n")
  blockMeshFile.write("    ("+str( c1)+" "+str( c1)+" "+str(height)+") //16\n")
  blockMeshFile.write("    ("+str(-c1)+" "+str( c1)+" "+str(height)+") //17\n")
  blockMeshFile.write("    ("+str(-c1)+" "+str(-c1)+" "+str(height)+") //18\n")
  blockMeshFile.write("    ("+str( c1)+" "+str(-c1)+" "+str(height)+") //19\n")
  blockMeshFile.write("    ("+str( c2)+" "+str( c2)+" "+str(height)+") //20\n")
  blockMeshFile.write("    ("+str(-c2)+" "+str( c2)+" "+str(height)+") //21\n")
  blockMeshFile.write("    ("+str(-c2)+" "+str(-c2)+" "+str(height)+") //22\n")
  blockMeshFile.write("    ("+str( c2)+" "+str(-c2)+" "+str(height)+") //23\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("blocks\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write("    hex ( 0  4  5  1 12 16 17 13) ("+str(nr1)+" "+str(nu)+" 1) simpleGrading ("+str(E1)+" 1 1)\n")
  blockMeshFile.write("    hex ( 1  5  6  2 13 17 18 14) ("+str(nr1)+" "+str(nu)+" 1) simpleGrading ("+str(E1)+" 1 1)\n")
  blockMeshFile.write("    hex ( 2  6  7  3 14 18 19 15) ("+str(nr1)+" "+str(nu)+" 1) simpleGrading ("+str(E1)+" 1 1)\n")
  blockMeshFile.write("    hex ( 3  7  4  0 15 19 16 12) ("+str(nr1)+" "+str(nu)+" 1) simpleGrading ("+str(E1)+" 1 1)\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    hex ( 4  8  9  5 16 20 21 17) ("+str(nr2)+" "+str(nu)+" 1) simpleGrading ("+str(E2)+" 1 1)\n")
  blockMeshFile.write("    hex ( 5  9 10  6 17 21 22 18) ("+str(nr2)+" "+str(nu)+" 1) simpleGrading ("+str(E2)+" 1 1)\n")
  blockMeshFile.write("    hex ( 6 10 11  7 18 22 23 19) ("+str(nr2)+" "+str(nu)+" 1) simpleGrading ("+str(E2)+" 1 1)\n")
  blockMeshFile.write("    hex ( 7 11  8  4 19 23 20 16) ("+str(nr2)+" "+str(nu)+" 1) simpleGrading ("+str(E2)+" 1 1)\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("edges\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write("    arc  0  1 (0 "+str(r0)+" 0)\n")
  blockMeshFile.write("    arc  1  2 ("+str(-r0)+" 0 0)\n")
  blockMeshFile.write("    arc  2  3 (0 "+str(-r0)+" 0)\n")
  blockMeshFile.write("    arc  3  0 ("+str(r0)+" 0 0)\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    arc  4  5 (0 "+str(r1)+" 0)\n")
  blockMeshFile.write("    arc  5  6 ("+str(-r1)+" 0 0)\n")
  blockMeshFile.write("    arc  6  7 (0 "+str(-r1)+" 0)\n")
  blockMeshFile.write("    arc  7  4 ("+str(r1)+" 0 0)\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    arc  8  9 (0 "+str(r2)+" 0)\n")
  blockMeshFile.write("    arc  9 10 ("+str(-r2)+" 0 0)\n")
  blockMeshFile.write("    arc 10 11 (0 "+str(-r2)+" 0)\n")
  blockMeshFile.write("    arc 11  8 ("+str(r2)+" 0 0)\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    arc 12 13 (0 "+str(r0)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 13 14 ("+str(-r0)+" 0 "+str(height)+")\n")
  blockMeshFile.write("    arc 14 15 (0 "+str(-r0)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 15 12 ("+str(r0)+" 0 "+str(height)+")\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    arc 16 17 (0 "+str(r1)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 17 18 ("+str(-r1)+" 0 "+str(height)+")\n")
  blockMeshFile.write("    arc 18 19 (0 "+str(-r1)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 19 16 ("+str(r1)+" 0 "+str(height)+")\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    arc 20 21 (0 "+str(r2)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 21 22 ("+str(-r2)+" 0 "+str(height)+")\n")
  blockMeshFile.write("    arc 22 23 (0 "+str(-r2)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 23 20 ("+str(r2)+" 0 "+str(height)+")\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("boundary\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write("    cylinder\n")
  blockMeshFile.write("    {\n")
  blockMeshFile.write("        type wall;\n")
  blockMeshFile.write("        faces\n")
  blockMeshFile.write("        (\n")
  blockMeshFile.write("            ( 0 12 13  1)\n")
  blockMeshFile.write("            ( 1 13 14  2)\n")
  blockMeshFile.write("            ( 2 14 15  3)\n")
  blockMeshFile.write("            ( 3 15 12  0)\n")
  blockMeshFile.write("        );\n")
  blockMeshFile.write("    }\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    inletOutlet")
  blockMeshFile.write("    {\n")
  blockMeshFile.write("        type patch;\n")
  blockMeshFile.write("        faces\n")
  blockMeshFile.write("        (\n")
  blockMeshFile.write("            ( 8  9 21 20)\n")
  blockMeshFile.write("            ( 9 10 22 21)\n")
  blockMeshFile.write("            (10 11 23 22)\n")
  blockMeshFile.write("            (11  8 20 23)\n")
  blockMeshFile.write("        );\n")
  blockMeshFile.write("    }\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    frontAndBack\n")
  blockMeshFile.write("    {\n")
  blockMeshFile.write("        type empty;\n")
  blockMeshFile.write("        faces\n")
  blockMeshFile.write("        (\n")
  blockMeshFile.write("            ( 0  1  5  4)\n")
  blockMeshFile.write("            ( 1  2  6  5)\n")
  blockMeshFile.write("            ( 2  3  7  6)\n")
  blockMeshFile.write("            ( 3  0  4  7)\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("            (16 17 13 12)\n")
  blockMeshFile.write("            (17 18 14 13)\n")
  blockMeshFile.write("            (18 19 15 14)\n")
  blockMeshFile.write("            (19 16 12 15)\n")
  blockMeshFile.write("        );\n")
  blockMeshFile.write("    }\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("mergePatchPairs\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("// ************************************************************************* //\n")
  blockMeshFile.close()
  runApp("blockMesh","",True,path)
  return
  
def createOversetBlockMesh(path, nu, nr, nbg, E, r0, r1, r2):
  """
  Create blockMesh with overset and background regions
  """
  if not os.path.exists(path+"/constant/polyMesh"):
    os.makedirs(path+"/constant/polyMesh")
  
  c0=math.sqrt(0.5)*r0
  c1=math.sqrt(0.5)*r1
  height=1.0

  blockMeshFile = open(path+"/constant/polyMesh/blockMeshDict",'w')
  blockMeshFile.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
  blockMeshFile.write("| =========                 |                                                 |\n")
  blockMeshFile.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
  blockMeshFile.write("|  \\    /   O peration     | Version:  2.1.x                                 |\n")
  blockMeshFile.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")
  blockMeshFile.write("|    \\/     M anipulation  |                                                 |\n")
  blockMeshFile.write("\*---------------------------------------------------------------------------*/\n")
  blockMeshFile.write("FoamFile\n")
  blockMeshFile.write("{\n")
  blockMeshFile.write("    version     2.0;\n")
  blockMeshFile.write("    format      ascii;\n")
  blockMeshFile.write("    class       dictionary;\n")
  blockMeshFile.write("    object      blockMeshDict;\n")
  blockMeshFile.write("}\n")
  blockMeshFile.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("convertToMeters 1;\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("vertices\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write("    ("+str( c0)+" "+str( c0)+" "+str( 0.0  )+") // 0\n")
  blockMeshFile.write("    ("+str(-c0)+" "+str( c0)+" "+str( 0.0  )+") // 1\n")
  blockMeshFile.write("    ("+str(-c0)+" "+str(-c0)+" "+str( 0.0  )+") // 2\n")
  blockMeshFile.write("    ("+str( c0)+" "+str(-c0)+" "+str( 0.0  )+") // 3\n")
  blockMeshFile.write("    ("+str( c1)+" "+str( c1)+" "+str( 0.0  )+") // 4\n")
  blockMeshFile.write("    ("+str(-c1)+" "+str( c1)+" "+str( 0.0  )+") // 5\n")
  blockMeshFile.write("    ("+str(-c1)+" "+str(-c1)+" "+str( 0.0  )+") // 6\n")
  blockMeshFile.write("    ("+str( c1)+" "+str(-c1)+" "+str( 0.0  )+") // 7\n")
  blockMeshFile.write("    ("+str(-r2)+" "+str(-r2)+" "+str( 0.0  )+") // 8\n")
  blockMeshFile.write("    ("+str( r2)+" "+str(-r2)+" "+str( 0.0  )+") // 9\n")
  blockMeshFile.write("    ("+str( r2)+" "+str( r2)+" "+str( 0.0  )+") //10\n")
  blockMeshFile.write("    ("+str(-r2)+" "+str( r2)+" "+str( 0.0  )+") //11\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    ("+str( c0)+" "+str( c0)+" "+str(height)+") //12\n")
  blockMeshFile.write("    ("+str(-c0)+" "+str( c0)+" "+str(height)+") //13\n")
  blockMeshFile.write("    ("+str(-c0)+" "+str(-c0)+" "+str(height)+") //14\n")
  blockMeshFile.write("    ("+str( c0)+" "+str(-c0)+" "+str(height)+") //15\n")
  blockMeshFile.write("    ("+str( c1)+" "+str( c1)+" "+str(height)+") //16\n")
  blockMeshFile.write("    ("+str(-c1)+" "+str( c1)+" "+str(height)+") //17\n")
  blockMeshFile.write("    ("+str(-c1)+" "+str(-c1)+" "+str(height)+") //18\n")
  blockMeshFile.write("    ("+str( c1)+" "+str(-c1)+" "+str(height)+") //19\n")
  blockMeshFile.write("    ("+str(-r2)+" "+str(-r2)+" "+str(height)+") //20\n")
  blockMeshFile.write("    ("+str( r2)+" "+str(-r2)+" "+str(height)+") //21\n")
  blockMeshFile.write("    ("+str( r2)+" "+str( r2)+" "+str(height)+") //22\n")
  blockMeshFile.write("    ("+str(-r2)+" "+str( r2)+" "+str(height)+") //23\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("blocks\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write("    hex ( 0  4  5  1 12 16 17 13) inner ("+str(nr)+" "+str(nu)+" 1) simpleGrading ("+str(E)+" 1 1)\n")
  blockMeshFile.write("    hex ( 1  5  6  2 13 17 18 14) inner ("+str(nr)+" "+str(nu)+" 1) simpleGrading ("+str(E)+" 1 1)\n")
  blockMeshFile.write("    hex ( 2  6  7  3 14 18 19 15) inner ("+str(nr)+" "+str(nu)+" 1) simpleGrading ("+str(E)+" 1 1)\n")
  blockMeshFile.write("    hex ( 3  7  4  0 15 19 16 12) inner ("+str(nr)+" "+str(nu)+" 1) simpleGrading ("+str(E)+" 1 1)\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    hex ( 8  9 10 11 20 21 22 23) outer ("+str(nbg)+" "+str(nbg)+" 1) simpleGrading (1 1 1)\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("edges\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write("    arc  0  1 (0 "+str(r0)+" 0)\n")
  blockMeshFile.write("    arc  1  2 ("+str(-r0)+" 0 0)\n")
  blockMeshFile.write("    arc  2  3 (0 "+str(-r0)+" 0)\n")
  blockMeshFile.write("    arc  3  0 ("+str(r0)+" 0 0)\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    arc  4  5 (0 "+str(r1)+" 0)\n")
  blockMeshFile.write("    arc  5  6 ("+str(-r1)+" 0 0)\n")
  blockMeshFile.write("    arc  6  7 (0 "+str(-r1)+" 0)\n")
  blockMeshFile.write("    arc  7  4 ("+str(r1)+" 0 0)\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    arc 12 13 (0 "+str(r0)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 13 14 ("+str(-r0)+" 0 "+str(height)+")\n")
  blockMeshFile.write("    arc 14 15 (0 "+str(-r0)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 15 12 ("+str(r0)+" 0 "+str(height)+")\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    arc 16 17 (0 "+str(r1)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 17 18 ("+str(-r1)+" 0 "+str(height)+")\n")
  blockMeshFile.write("    arc 18 19 (0 "+str(-r1)+" "+str(height)+")\n")
  blockMeshFile.write("    arc 19 16 ("+str(r1)+" 0 "+str(height)+")\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("boundary\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write("    cylinder\n")
  blockMeshFile.write("    {\n")
  blockMeshFile.write("        type wall;\n")
  blockMeshFile.write("        faces\n")
  blockMeshFile.write("        (\n")
  blockMeshFile.write("            ( 0 12 13  1)\n")
  blockMeshFile.write("            ( 1 13 14  2)\n")
  blockMeshFile.write("            ( 2 14 15  3)\n")
  blockMeshFile.write("            ( 3 15 12  0)\n")
  blockMeshFile.write("        );\n")
  blockMeshFile.write("    }\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    inner")
  blockMeshFile.write("    {\n")
  blockMeshFile.write("        type bellerophon;\n")
  blockMeshFile.write("        donorZone outer;\n")
  blockMeshFile.write("        faces\n")
  blockMeshFile.write("        (\n")
  blockMeshFile.write("            ( 4  5 17 16)\n")
  blockMeshFile.write("            ( 5  6 18 17)\n")
  blockMeshFile.write("            ( 6  7 19 18)\n")
  blockMeshFile.write("            ( 7  4 16 19)\n")
  blockMeshFile.write("        );\n")
  blockMeshFile.write("    }\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    inlet")
  blockMeshFile.write("    {\n")
  blockMeshFile.write("        type patch;\n")
  blockMeshFile.write("        faces\n")
  blockMeshFile.write("        (\n")
  blockMeshFile.write("            ( 8 20 23 11)\n")
  blockMeshFile.write("        );\n")
  blockMeshFile.write("    }\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    outlet")
  blockMeshFile.write("    {\n")
  blockMeshFile.write("        type patch;\n")
  blockMeshFile.write("        faces\n")
  blockMeshFile.write("        (\n")
  blockMeshFile.write("            ( 9 10 22 21)\n")
  blockMeshFile.write("        );\n")
  blockMeshFile.write("    }\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    sides")
  blockMeshFile.write("    {\n")
  blockMeshFile.write("        type patch;\n")
  blockMeshFile.write("        faces\n")
  blockMeshFile.write("        (\n")
  blockMeshFile.write("            (10 11 23 22)\n")
  blockMeshFile.write("            ( 8  9 21 20)\n")
  blockMeshFile.write("        );\n")
  blockMeshFile.write("    }\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("    frontAndBack\n")
  blockMeshFile.write("    {\n")
  blockMeshFile.write("        type empty;\n")
  blockMeshFile.write("        faces\n")
  blockMeshFile.write("        (\n")
  blockMeshFile.write("            ( 0  1  5  4)\n")
  blockMeshFile.write("            ( 1  2  6  5)\n")
  blockMeshFile.write("            ( 2  3  7  6)\n")
  blockMeshFile.write("            ( 3  0  4  7)\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("            ( 8  9 10 11)\n")
  blockMeshFile.write("            (23 22 21 20)\n")
  blockMeshFile.write("        );\n")
  blockMeshFile.write("    }\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("mergePatchPairs\n")
  blockMeshFile.write("(\n")
  blockMeshFile.write(");\n")
  blockMeshFile.write("\n")
  blockMeshFile.write("// ************************************************************************* //\n")
  blockMeshFile.close()
  runApp("blockMesh","",True,path)
  return
  
def runApp(appName, flags="", log=True, path = None):
  """
  run application with flags/options and write output to log.appName if log is
  enabled
  """
  if not path:
    path = os.getcwd()
    oldPath = path
  else:
    oldPath = os.getcwd()
  os.chdir(path)
  print("  Running "+appName)
  if log:
    renameLog("log."+appName)
    os.system(appName+" "+flags+" 2>&1 > log."+appName)
  else:
    os.system(appName+" "+flags)
  os.chdir(oldPath)
  return

def runParallel(appName, nProcs, flags="", log=True, path = None):
  """
  run application with flags/options and write output to log.appName if log is
  enabled
  """
  if not path:
    path = os.getcwd()
    oldPath = path
  else:
    oldPath = os.getcwd()
  os.chdir(path)
  print("  Running "+appName)
  if log:
    renameLog("log."+appName)
    os.system("mpirun -np "+str(nProcs)+" "+appName+" "+flags+" 2>&1 > log."+appName)
  else:
    os.system("mpirun -np "+str(nProcs)+" "+appName+" "+flags)
  os.chdir(oldPath)
  return

def decompAndRun(appName, nProcs, args="", onlyLatest=False):
  """
  Decompose, run solver and reconstruct
  """
  runApp("decomposePar")
  runParallel(appName,nProcs,"-parallel "+args)
  if onlyLatest:
    runApp("reconstructPar","-latestTime")
  else:
    runApp("reconstructPar")    

def renameLog(fileName):
  """
  move log file
  """
  if os.path.isfile(fileName):
    i=0
    while os.path.isfile(fileName+"."+str(i)):
      i=i+1
  return
      
def extractFluxConservation(caseName, logName, resultName):
  """
  Extract flux conservation error from caseName/logFile and write it to
  resultName
  """
  logFile = open(caseName+"/"+logName,'r')
  resultFile = open(resultName,'w')
  resultFile.write("#time\terror_0\terror_1\t...")
  for line in logFile:
     if line.startswith("Time = "):
       t = line.split()[2]
       resultFile.write("\n"+t)
     elif line.startswith("Interface error : "):
       lineArray=line.split()
       if len(lineArray) is not 6:
         print("Error line doesn't match: "+line)
         raise ValueError
       else:
         resultFile.write("\t"+str(float(lineArray[3])/float(lineArray[5])))
  logFile.close()
  resultFile.close()
