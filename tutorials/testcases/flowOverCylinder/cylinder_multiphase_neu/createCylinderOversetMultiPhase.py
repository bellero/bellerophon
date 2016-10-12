#!/usr/bin/env python3

##modules
import os
import math
pi=math.pi;
cos=math.cos;

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

##geometry 
D=1.0;           #diameter cylinder 
Do=2.0;          #overset domain diameter

T=4.0;           #water depth
H=2.0;           #height above sea level

xup=6.0;         #width domain in x-direction (in front of the cylinder)
xdown=15.0;      #width domain in x-direction (behind the cylinder)
width=12.0;      #width domain in y-direction (left side) 

hinterface=0.4;  #interface refinement zone in z-direction
ainterface=5.0   # aspect ratio at interface (cell size at outer overset interface)

##meshing
#cells
nc=50  	         #cells in circumferential direction
aspect=2.5       #aspect ratio of cells at cylinder surface

#grading
epso=1.07;	      #expansion ratio in radial direction

epsbg=1.1;	      #expansion ratio in horizontal directions in background

epsz=1.1;	      #expansion ratio in vertical direction in background

print("********** write blockMeshDict **********")

xmin = -xup * D
xmax = xdown * D
ymax = 0.5 * width * D
ymin = -ymax
zimax = 0.5 * hinterface * D
zimin = -zimax
zmin = -T * D
zmax = H * D

r = 0.5 * D
ro = 0.5 * Do
deltaU0 = 0.5 * pi * r / (aspect * nc)  # Height of cells on cylinder surface
deltaU1 = 0.5 * pi * ro / nc            # Height of cells on block boundary
deltaZ = deltaU1 / ainterface

c0 = math.sqrt(0.5)*r
c1 = math.sqrt(0.5)*ro

Eo = deltaU1 / deltaU0            # Grading of overset domain
nr = nFromLEDelta(ro-r, Eo, deltaU0) # Number of cells in radial direction
nz = int(abs(zimax-zimin) / deltaZ)

rinner = ro + 5.0 * deltaU1       # Size of inner block in background
ninner = int(2.0 * rinner / deltaU1)    # cells in inner block in background

nup = nFromLEpsilonDelta(abs(abs(xmin)-rinner), epsbg, deltaU1*epsbg)
ndown = nFromLEpsilonDelta(abs(abs(xmax)-rinner), epsbg, deltaU1*epsbg)
nside = nFromLEpsilonDelta(abs(abs(ymin)-rinner), epsbg, deltaU1*epsbg)
ntop = nFromLEpsilonDelta(abs(abs(zmax)-abs(zimax)), epsz, deltaZ*epsz)
nbottom = nFromLEpsilonDelta(abs(abs(zmin)-abs(zimin)), epsz, deltaZ*epsz)

Eup = 1.0 / (epsbg**(nup-1))
Edown = epsbg**(ndown-1)
Eside = epsbg**(nside-1)
Etop = epsz**(ntop-1)
Ebottom = epsz**(nbottom-1)

if not os.path.exists("3d/constant/polyMesh"):
  os.makedirs("3d/constant/polyMesh")

blockMeshFile = open("3d/constant/polyMesh/blockMeshDict",'w')
blockMeshFile.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
blockMeshFile.write("| =========                 |                                                 |\n")
blockMeshFile.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
blockMeshFile.write("|  \\    /   O peration     | Version:  2.1.x                                 |\n")
blockMeshFile.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")
blockMeshFile.write("|    \\/     M anipulation  |                                                 |\n")
blockMeshFile.write("\\*---------------------------------------------------------------------------*/\n")
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
blockMeshFile.write("    ("+str( c0)+" "+str( c0)+" "+str(zimin)+") // 0\n")
blockMeshFile.write("    ("+str(-c0)+" "+str( c0)+" "+str(zimin)+") // 1\n")
blockMeshFile.write("    ("+str(-c0)+" "+str(-c0)+" "+str(zimin)+") // 2\n")
blockMeshFile.write("    ("+str( c0)+" "+str(-c0)+" "+str(zimin)+") // 3\n")
blockMeshFile.write("    ("+str( c1)+" "+str( c1)+" "+str(zimin)+") // 4\n")
blockMeshFile.write("    ("+str(-c1)+" "+str( c1)+" "+str(zimin)+") // 5\n")
blockMeshFile.write("    ("+str(-c1)+" "+str(-c1)+" "+str(zimin)+") // 6\n")
blockMeshFile.write("    ("+str( c1)+" "+str(-c1)+" "+str(zimin)+") // 7\n")
blockMeshFile.write("\n")
blockMeshFile.write("    ("+str( c0)+" "+str( c0)+" "+str(zimax)+") // 8\n")
blockMeshFile.write("    ("+str(-c0)+" "+str( c0)+" "+str(zimax)+") // 9\n")
blockMeshFile.write("    ("+str(-c0)+" "+str(-c0)+" "+str(zimax)+") //10\n")
blockMeshFile.write("    ("+str( c0)+" "+str(-c0)+" "+str(zimax)+") //11\n")
blockMeshFile.write("    ("+str( c1)+" "+str( c1)+" "+str(zimax)+") //12\n")
blockMeshFile.write("    ("+str(-c1)+" "+str( c1)+" "+str(zimax)+") //13\n")
blockMeshFile.write("    ("+str(-c1)+" "+str(-c1)+" "+str(zimax)+") //14\n")
blockMeshFile.write("    ("+str( c1)+" "+str(-c1)+" "+str(zimax)+") //15\n")
blockMeshFile.write("\n")
blockMeshFile.write("    ("+str(  xmin )+" "+str(  ymax )+" "+str(zimin)+") //16\n")
blockMeshFile.write("    ("+str(  xmin )+" "+str( rinner)+" "+str(zimin)+") //17\n")
blockMeshFile.write("    ("+str(  xmin )+" "+str(-rinner)+" "+str(zimin)+") //18\n")
blockMeshFile.write("    ("+str(  xmin )+" "+str(  ymin )+" "+str(zimin)+") //19\n")
blockMeshFile.write("    ("+str(-rinner)+" "+str(  ymax )+" "+str(zimin)+") //20\n")
blockMeshFile.write("    ("+str(-rinner)+" "+str( rinner)+" "+str(zimin)+") //21\n")
blockMeshFile.write("    ("+str(-rinner)+" "+str(-rinner)+" "+str(zimin)+") //22\n")
blockMeshFile.write("    ("+str(-rinner)+" "+str(  ymin )+" "+str(zimin)+") //23\n")
blockMeshFile.write("    ("+str( rinner)+" "+str(  ymax )+" "+str(zimin)+") //24\n")
blockMeshFile.write("    ("+str( rinner)+" "+str( rinner)+" "+str(zimin)+") //25\n")
blockMeshFile.write("    ("+str( rinner)+" "+str(-rinner)+" "+str(zimin)+") //26\n")
blockMeshFile.write("    ("+str( rinner)+" "+str(  ymin )+" "+str(zimin)+") //27\n")
blockMeshFile.write("    ("+str(  xmax )+" "+str(  ymax )+" "+str(zimin)+") //28\n")
blockMeshFile.write("    ("+str(  xmax )+" "+str( rinner)+" "+str(zimin)+") //29\n")
blockMeshFile.write("    ("+str(  xmax )+" "+str(-rinner)+" "+str(zimin)+") //30\n")
blockMeshFile.write("    ("+str(  xmax )+" "+str(  ymin )+" "+str(zimin)+") //31\n")
blockMeshFile.write("\n")
blockMeshFile.write("    ("+str(  xmin )+" "+str(  ymax )+" "+str(zimax)+") //32\n")
blockMeshFile.write("    ("+str(  xmin )+" "+str( rinner)+" "+str(zimax)+") //33\n")
blockMeshFile.write("    ("+str(  xmin )+" "+str(-rinner)+" "+str(zimax)+") //34\n")
blockMeshFile.write("    ("+str(  xmin )+" "+str(  ymin )+" "+str(zimax)+") //35\n")
blockMeshFile.write("    ("+str(-rinner)+" "+str(  ymax )+" "+str(zimax)+") //36\n")
blockMeshFile.write("    ("+str(-rinner)+" "+str( rinner)+" "+str(zimax)+") //37\n")
blockMeshFile.write("    ("+str(-rinner)+" "+str(-rinner)+" "+str(zimax)+") //38\n")
blockMeshFile.write("    ("+str(-rinner)+" "+str(  ymin )+" "+str(zimax)+") //39\n")
blockMeshFile.write("    ("+str( rinner)+" "+str(  ymax )+" "+str(zimax)+") //40\n")
blockMeshFile.write("    ("+str( rinner)+" "+str( rinner)+" "+str(zimax)+") //41\n")
blockMeshFile.write("    ("+str( rinner)+" "+str(-rinner)+" "+str(zimax)+") //42\n")
blockMeshFile.write("    ("+str( rinner)+" "+str(  ymin )+" "+str(zimax)+") //43\n")
blockMeshFile.write("    ("+str(  xmax )+" "+str(  ymax )+" "+str(zimax)+") //44\n")
blockMeshFile.write("    ("+str(  xmax )+" "+str( rinner)+" "+str(zimax)+") //45\n")
blockMeshFile.write("    ("+str(  xmax )+" "+str(-rinner)+" "+str(zimax)+") //46\n")
blockMeshFile.write("    ("+str(  xmax )+" "+str(  ymin )+" "+str(zimax)+") //47\n")
blockMeshFile.write("\n")
blockMeshFile.write(");\n")
blockMeshFile.write("\n")
blockMeshFile.write("blocks\n")
blockMeshFile.write("(\n")
blockMeshFile.write("    hex ( 0  4  5  1  8 12 13  9) inner ("+str(nr)+" "+str(nc)+" "+str(nz)+") simpleGrading ("+str(Eo)+" 1 1)\n")
blockMeshFile.write("    hex ( 1  5  6  2  9 13 14 10) inner ("+str(nr)+" "+str(nc)+" "+str(nz)+") simpleGrading ("+str(Eo)+" 1 1)\n")
blockMeshFile.write("    hex ( 2  6  7  3 10 14 15 11) inner ("+str(nr)+" "+str(nc)+" "+str(nz)+") simpleGrading ("+str(Eo)+" 1 1)\n")
blockMeshFile.write("    hex ( 3  7  4  0 11 15 12  8) inner ("+str(nr)+" "+str(nc)+" "+str(nz)+") simpleGrading ("+str(Eo)+" 1 1)\n")
blockMeshFile.write("\n")
blockMeshFile.write("    hex (17 21 20 16 33 37 36 32) outer ("+str(  nup )+" "+str( nside)+" "+str(nz)+") simpleGrading ("+str(  Eup  )+" "+str(  Eside  )+" 1)\n")
blockMeshFile.write("    hex (18 22 21 17 34 38 37 33) outer ("+str(  nup )+" "+str(ninner)+" "+str(nz)+") simpleGrading ("+str(  Eup  )+" "+str(    1    )+" 1)\n")
blockMeshFile.write("    hex (19 23 22 18 35 39 38 34) outer ("+str(  nup )+" "+str( nside)+" "+str(nz)+") simpleGrading ("+str(  Eup  )+" "+str(1.0/Eside)+" 1)\n")
blockMeshFile.write("\n")
blockMeshFile.write("    hex (21 25 24 20 37 41 40 36) outer ("+str(ninner)+" "+str( nside)+" "+str(nz)+") simpleGrading ("+str(   1   )+" "+str(  Eside  )+" 1)\n")
blockMeshFile.write("    hex (22 26 25 21 38 42 41 37) outer ("+str(ninner)+" "+str(ninner)+" "+str(nz)+") simpleGrading ("+str(   1   )+" "+str(    1    )+" 1)\n")
blockMeshFile.write("    hex (23 27 26 22 39 43 42 38) outer ("+str(ninner)+" "+str( nside)+" "+str(nz)+") simpleGrading ("+str(   1   )+" "+str(1.0/Eside)+" 1)\n")
blockMeshFile.write("\n")
blockMeshFile.write("    hex (25 29 28 24 41 45 44 40) outer ("+str( ndown)+" "+str( nside)+" "+str(nz)+") simpleGrading ("+str( Edown )+" "+str(  Eside  )+" 1)\n")
blockMeshFile.write("    hex (26 30 29 25 42 46 45 41) outer ("+str( ndown)+" "+str(ninner)+" "+str(nz)+") simpleGrading ("+str( Edown )+" "+str(    1    )+" 1)\n")
blockMeshFile.write("    hex (27 31 30 26 43 47 46 42) outer ("+str( ndown)+" "+str( nside)+" "+str(nz)+") simpleGrading ("+str( Edown )+" "+str(1.0/Eside)+" 1)\n")
blockMeshFile.write(");\n")
blockMeshFile.write("\n")
blockMeshFile.write("edges\n")
blockMeshFile.write("(\n")
blockMeshFile.write("    arc  0  1 (0 "+str(r)+" "+str(zimin)+")\n")
blockMeshFile.write("    arc  1  2 ("+str(-r)+" 0 "+str(zimin)+")\n")
blockMeshFile.write("    arc  2  3 (0 "+str(-r)+" "+str(zimin)+")\n")
blockMeshFile.write("    arc  3  0 ("+str(r)+" 0 "+str(zimin)+")\n")
blockMeshFile.write("\n")
blockMeshFile.write("    arc  4  5 (0 "+str(ro)+" "+str(zimin)+")\n")
blockMeshFile.write("    arc  5  6 ("+str(-ro)+" 0 "+str(zimin)+")\n")
blockMeshFile.write("    arc  6  7 (0 "+str(-ro)+" "+str(zimin)+")\n")
blockMeshFile.write("    arc  7  4 ("+str(ro)+" 0 "+str(zimin)+")\n")
blockMeshFile.write("\n")
blockMeshFile.write("    arc  8  9 (0 "+str(r)+" "+str(zimax)+")\n")
blockMeshFile.write("    arc  9 10 ("+str(-r)+" 0 "+str(zimax)+")\n")
blockMeshFile.write("    arc 10 11 (0 "+str(-r)+" "+str(zimax)+")\n")
blockMeshFile.write("    arc 11  8 ("+str(r)+" 0 "+str(zimax)+")\n")
blockMeshFile.write("\n")
blockMeshFile.write("    arc 12 13 (0 "+str(ro)+" "+str(zimax)+")\n")
blockMeshFile.write("    arc 13 14 ("+str(-ro)+" 0 "+str(zimax)+")\n")
blockMeshFile.write("    arc 14 15 (0 "+str(-ro)+" "+str(zimax)+")\n")
blockMeshFile.write("    arc 15 12 ("+str(ro)+" 0 "+str(zimax)+")\n")
blockMeshFile.write(");\n")
blockMeshFile.write("\n")
blockMeshFile.write("boundary\n")
blockMeshFile.write("(\n")
blockMeshFile.write("    cylinder\n")
blockMeshFile.write("    {\n")
blockMeshFile.write("        type wall;\n")
blockMeshFile.write("        faces\n")
blockMeshFile.write("        (\n")
blockMeshFile.write("            ( 0  8  9  1)\n")
blockMeshFile.write("            ( 1  9 10  2)\n")
blockMeshFile.write("            ( 2 10 11  3)\n")
blockMeshFile.write("            ( 3 11  8  0)\n")
blockMeshFile.write("        );\n")
blockMeshFile.write("    }\n")
blockMeshFile.write("\n")
blockMeshFile.write("    inner\n")
blockMeshFile.write("    {\n")
blockMeshFile.write("        type bellerophon;\n")
blockMeshFile.write("        donorZone outer;\n")
blockMeshFile.write("        faces\n")
blockMeshFile.write("        (\n")
blockMeshFile.write("            ( 4  5 13 12)\n")
blockMeshFile.write("            ( 5  6 14 13)\n")
blockMeshFile.write("            ( 6  7 15 14)\n")
blockMeshFile.write("            ( 7  4 12 15)\n")
blockMeshFile.write("        );\n")
blockMeshFile.write("    }\n")
blockMeshFile.write("\n")
blockMeshFile.write("    inlet\n")
blockMeshFile.write("    {\n")
blockMeshFile.write("        type patch;\n")
blockMeshFile.write("        faces\n")
blockMeshFile.write("        (\n")
blockMeshFile.write("            (16 17 33 32)\n")
blockMeshFile.write("            (17 18 34 33)\n")
blockMeshFile.write("            (18 19 35 34)\n")
blockMeshFile.write("        );\n")
blockMeshFile.write("    }\n")
blockMeshFile.write("\n")
blockMeshFile.write("    outlet\n")
blockMeshFile.write("    {\n")
blockMeshFile.write("        type patch;\n")
blockMeshFile.write("        faces\n")
blockMeshFile.write("        (\n")
blockMeshFile.write("            (29 28 44 45)\n")
blockMeshFile.write("            (30 29 45 46)\n")
blockMeshFile.write("            (31 30 46 47)\n")
blockMeshFile.write("        );\n")
blockMeshFile.write("    }\n")
blockMeshFile.write("\n")
blockMeshFile.write("    sides\n")
blockMeshFile.write("    {\n")
blockMeshFile.write("        type patch;\n")
blockMeshFile.write("        faces\n")
blockMeshFile.write("        (\n")
blockMeshFile.write("            (20 16 32 36)\n")
blockMeshFile.write("            (24 20 36 40)\n")
blockMeshFile.write("            (28 24 40 44)\n")
blockMeshFile.write("\n")
blockMeshFile.write("            (19 23 39 35)\n")
blockMeshFile.write("            (23 27 43 39)\n")
blockMeshFile.write("            (27 31 47 43)\n")
blockMeshFile.write("        );\n")
blockMeshFile.write("    }\n")
blockMeshFile.write("\n")
blockMeshFile.write("    top\n")
blockMeshFile.write("    {\n")
blockMeshFile.write("        type patch;\n")
blockMeshFile.write("        faces\n")
blockMeshFile.write("        (\n")
blockMeshFile.write("            ( 12 13  9  8)\n")
blockMeshFile.write("            ( 13 14 10  9)\n")
blockMeshFile.write("            ( 14 15 11 10)\n")
blockMeshFile.write("            ( 15 12  8 11)\n")
blockMeshFile.write("\n")
blockMeshFile.write("            (33 37 36 32)\n")
blockMeshFile.write("            (34 38 37 33)\n")
blockMeshFile.write("            (35 39 38 34)\n")
blockMeshFile.write("\n")
blockMeshFile.write("            (37 41 40 36)\n")
blockMeshFile.write("            (38 42 41 37)\n")
blockMeshFile.write("            (39 43 42 38)\n")
blockMeshFile.write("\n")
blockMeshFile.write("            (41 45 44 40)\n")
blockMeshFile.write("            (42 46 45 41)\n")
blockMeshFile.write("            (43 47 46 42)\n")
blockMeshFile.write("        );\n")
blockMeshFile.write("    }\n")
blockMeshFile.write("\n")
blockMeshFile.write("    bottom\n")
blockMeshFile.write("    {\n")
blockMeshFile.write("        type patch;\n")
blockMeshFile.write("        faces\n")
blockMeshFile.write("        (\n")
blockMeshFile.write("            ( 0  1  5  4)\n")
blockMeshFile.write("            ( 1  2  6  5)\n")
blockMeshFile.write("            ( 2  3  7  6)\n")
blockMeshFile.write("            ( 3  0  4  7)\n")
blockMeshFile.write("\n")
blockMeshFile.write("            (16 20 21 17)\n")
blockMeshFile.write("            (17 21 22 18)\n")
blockMeshFile.write("            (18 22 23 19)\n")
blockMeshFile.write("\n")
blockMeshFile.write("            (20 24 25 21)\n")
blockMeshFile.write("            (21 25 26 22)\n")
blockMeshFile.write("            (22 26 27 23)\n")
blockMeshFile.write("\n")
blockMeshFile.write("            (24 28 29 25)\n")
blockMeshFile.write("            (25 29 30 26)\n")
blockMeshFile.write("            (26 30 31 27)\n")
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

print("********* write extrudeMeshDict **********")
extrudeMeshDict=open("3d/system/extrudeMeshDict.bottom", "w")
extrudeMeshDict.write("/*--------------------------------*- C++ -*----------------------------------*\\ \n");
extrudeMeshDict.write("| =========                |                                                  | \n");
extrudeMeshDict.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n");
extrudeMeshDict.write("|  \\    /   O peration     | Version:  2.3.0                                 | \n");
extrudeMeshDict.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n");
extrudeMeshDict.write("|    \\/     M anipulation  |                                                 | \n");
extrudeMeshDict.write("\\*---------------------------------------------------------------------------*/ \n");
extrudeMeshDict.write("FoamFile\n");
extrudeMeshDict.write("{\n");
extrudeMeshDict.write("    version     2.0;\n");
extrudeMeshDict.write("    format      ascii;\n");
extrudeMeshDict.write("    class       dictionary;\n");
extrudeMeshDict.write("    object      extrudeProperties;\n");
extrudeMeshDict.write("}\n");
extrudeMeshDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n");
extrudeMeshDict.write(" \n");
extrudeMeshDict.write("constructFrom mesh;\n");
extrudeMeshDict.write('sourceCase ".";\n');
extrudeMeshDict.write("sourcePatches (bottom);\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("flipNormals false;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("extrudeModel        linearDirection;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("nLayers             "+str(nbottom)+";\n");
extrudeMeshDict.write("expansionRatio      "+str(epsz)+";\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("linearDirectionCoeffs\n");
extrudeMeshDict.write("{ \n");
extrudeMeshDict.write("    axisPt          (0 0 0); \n");
extrudeMeshDict.write("    direction       (0 0 -1); \n");
extrudeMeshDict.write("    thickness       "+str(abs(zmin-zimin))+";\n");
extrudeMeshDict.write("}\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("mergeFaces false;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("mergeTol 0;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n");
extrudeMeshDict.close()

extrudeMeshDict=open("3d/system/extrudeMeshDict.top", "w")
extrudeMeshDict.write("/*--------------------------------*- C++ -*----------------------------------*\\ \n");
extrudeMeshDict.write("| =========                |                                                  | \n");
extrudeMeshDict.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n");
extrudeMeshDict.write("|  \\    /   O peration     | Version:  2.3.0                                 | \n");
extrudeMeshDict.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n");
extrudeMeshDict.write("|    \\/     M anipulation  |                                                 | \n");
extrudeMeshDict.write("\\*---------------------------------------------------------------------------*/ \n");
extrudeMeshDict.write("FoamFile\n");
extrudeMeshDict.write("{\n");
extrudeMeshDict.write("    version     2.0;\n");
extrudeMeshDict.write("    format      ascii;\n");
extrudeMeshDict.write("    class       dictionary;\n");
extrudeMeshDict.write("    object      extrudeProperties;\n");
extrudeMeshDict.write("}\n");
extrudeMeshDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n");
extrudeMeshDict.write(" \n");
extrudeMeshDict.write("constructFrom mesh;\n");
extrudeMeshDict.write('sourceCase ".";\n');
extrudeMeshDict.write("sourcePatches (top);\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("flipNormals false;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("extrudeModel        linearDirection;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("nLayers             "+str(ntop)+";\n");
extrudeMeshDict.write("expansionRatio      "+str(epsz)+";\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("linearDirectionCoeffs\n");
extrudeMeshDict.write("{ \n");
extrudeMeshDict.write("    axisPt          (0 0 0); \n");
extrudeMeshDict.write("    direction       (0 0 1); \n");
extrudeMeshDict.write("    thickness       "+str(abs(zmax-zimax))+";\n");
extrudeMeshDict.write("}\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("mergeFaces false;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("mergeTol 0;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n");
extrudeMeshDict.close()

extrudeMeshDict=open("2d/system/extrudeMeshDict", "w")
extrudeMeshDict.write("/*--------------------------------*- C++ -*----------------------------------*\\ \n");
extrudeMeshDict.write("| =========                 |                                                 | \n");
extrudeMeshDict.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           | \n");
extrudeMeshDict.write("|  \\    /   O peration     | Version:  2.3.0                                 | \n");
extrudeMeshDict.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      | \n");
extrudeMeshDict.write("|    \\/     M anipulation  |                                                 | \n");
extrudeMeshDict.write("\\*---------------------------------------------------------------------------*/ \n");
extrudeMeshDict.write("FoamFile\n");
extrudeMeshDict.write("{\n");
extrudeMeshDict.write("    version     2.0;\n");
extrudeMeshDict.write("    format      ascii;\n");
extrudeMeshDict.write("    class       dictionary;\n");
extrudeMeshDict.write("    object      extrudeProperties;\n");
extrudeMeshDict.write("}\n");
extrudeMeshDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // \n");
extrudeMeshDict.write(" \n");
extrudeMeshDict.write("constructFrom patch;\n");
extrudeMeshDict.write('sourceCase "../3d";\n');
extrudeMeshDict.write("sourcePatches ( bottom );\n");
extrudeMeshDict.write("exposedPatchName top;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("flipNormals true;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("extrudeModel        linearDirection;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("nLayers             1;\n");
extrudeMeshDict.write("expansionRatio      1;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("linearDirectionCoeffs\n");
extrudeMeshDict.write("{ \n");
extrudeMeshDict.write("    axisPt          (0 0 0); \n");
extrudeMeshDict.write("    direction       (0 0 1); \n");
extrudeMeshDict.write("    thickness       "+str(abs(zmax-zmin))+";\n");
extrudeMeshDict.write("}\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("mergeFaces false;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("mergeTol 0;\n");
extrudeMeshDict.write("\n");
extrudeMeshDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n");
extrudeMeshDict.close()

topoSetDict = open("3d/system/topoSetDict","w")
topoSetDict.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
topoSetDict.write("| =========                 |                                                 |\n")   
topoSetDict.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")   
topoSetDict.write("|  \\    /   O peration     | Version:  2.3.0                                 |\n")   
topoSetDict.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")   
topoSetDict.write("|    \\/     M anipulation  |                                                 |\n")   
topoSetDict.write("\\*---------------------------------------------------------------------------*/\n")
topoSetDict.write("FoamFile\n")
topoSetDict.write("{\n")
topoSetDict.write("    version     2.0;\n")
topoSetDict.write("    format      ascii;\n")
topoSetDict.write("    class       dictionary;\n")
topoSetDict.write("    object      topoSetDict;\n")
topoSetDict.write("}\n")
topoSetDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
topoSetDict.write("\n")
topoSetDict.write("actions\n")
topoSetDict.write("(\n")
topoSetDict.write("    {\n")
topoSetDict.write("        name    holeSource;\n")
topoSetDict.write("        type    cellSet;\n")
topoSetDict.write("        action  new;\n")
topoSetDict.write("        source  cylinderToCell;\n")
topoSetDict.write("        sourceInfo\n")
topoSetDict.write("        {\n")
topoSetDict.write("            p1      (0 0 "+str(zmin)+");\n")
topoSetDict.write("            p2      (0 0 "+str(zmax)+");\n")
topoSetDict.write("            radius  "+str(ro)+";\n")
topoSetDict.write("        }\n")
topoSetDict.write("    }\n")
topoSetDict.write(");\n")
topoSetDict.write("\n")
topoSetDict.write("// ************************************************************************* //\n")
topoSetDict.close()
