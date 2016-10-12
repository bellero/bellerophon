# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 19:51:16 2015

@author: andreas
"""

import numpy as np
import os

def makeMesh(z=4):
    Ro=1.0   # Height of Domain
    L=1.6   # Length of Domain
    
    Rb=0.5   # Radius of Blades
    Lb=0.2   # Length on Blades
    Rs=0.2   # Radius of Shaft
    Do=0.3   # Distance for Overset Region
    
    Nou=10   # Number of cells for Overset
    
    Du=Do/float(Nou)      # Cell Size at Interface
    
    print("Cell size at interface: "+str(Du))
    
    Nbx=int(L/Du/1.1)
    Nbr=int((Ro-Rs)/Du/1.1)
    Nbt=int(0.5*np.pi*(Rb+Do)/Du/1.1)
    
    Noxo=int(Do/Du)
    Noxi=int(Lb/Du)
    Nori=int(Rb/Du)
    Noro=int(Do/Du)
    Not=int(2.0*np.pi*(Rb+Do)/Du/z)
    
    os.system("rm -rf constant/polyMesh 0")
    os.mkdir("0")
    
    blockMeshDict = open("system/blockMeshDict","w")
    blockMeshDict.write("/*--------------------------------*- C++ -*----------------------------------*\\\n");
    blockMeshDict.write("| =========                 |                                                 |\n");
    blockMeshDict.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n");
    blockMeshDict.write("|  \\\\    /   O peration     | Version:  dev                                   |\n");
    blockMeshDict.write("|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n");
    blockMeshDict.write("|    \\\\/     M anipulation  |                                                 |\n");
    blockMeshDict.write("\\*---------------------------------------------------------------------------*/\n");
    blockMeshDict.write("FoamFile\n");
    blockMeshDict.write("{\n");
    blockMeshDict.write("    version     2.0;\n");
    blockMeshDict.write("    format      ascii;\n");
    blockMeshDict.write("    class       dictionary;\n");
    blockMeshDict.write("    object      blockMeshDict;\n");
    blockMeshDict.write("}\n");
    blockMeshDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n");
    blockMeshDict.write("\n");
    blockMeshDict.write("convertToMeters 1;\n");
    blockMeshDict.write("\n");
    blockMeshDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n");
    blockMeshDict.write("\n");
    blockMeshDict.write("vertices\n");
    blockMeshDict.write("(\n");
    for i in range(4):
        alpha=0.5*np.pi*i
        blockMeshDict.write("    ("+str(-0.5*L)+" "+str(Rs*np.sin(alpha))+" "+str(Rs*np.cos(alpha))+")\n");
        blockMeshDict.write("    ("+str(0.5*L)+" "+str(Rs*np.sin(alpha))+" "+str(Rs*np.cos(alpha))+")\n");
        blockMeshDict.write("    ("+str(0.5*L)+" "+str(Ro*np.sin(alpha))+" "+str(Ro*np.cos(alpha))+")\n");
        blockMeshDict.write("    ("+str(-0.5*L)+" "+str(Ro*np.sin(alpha))+" "+str(Ro*np.cos(alpha))+")\n");
    for i in range(z):
      alpha=2.0*np.pi/z*i
      for rcoord in [Rs, Rb,Rb+Do]:
        for x in [-Do-0.5*Lb, -0.5*Lb, 0.5*Lb, Do+0.5*Lb]:
          xcoord=str(x)
          blockMeshDict.write("    ("+xcoord+" "+str(rcoord*np.sin(alpha))+" "+str(rcoord*np.cos(alpha))+")\n");
    blockMeshDict.write(");\n");
    blockMeshDict.write("\n");
    blockMeshDict.write("blocks\n");
    blockMeshDict.write("(\n");
    blockMeshDict.write("    // block0\n");
    for i in range(4):
        blockMeshDict.write("    hex ("+" ".join([str(int(x+4*i)%16) for x in [0,1,5,4,3,2,6,7]])+") stator ("+str(Nbx)+" "+str(Nbt)+" "+str(Nbr)+") simpleGrading (1 1 1)\n");
    for i in range(z):
      for j in range(3):
        for k in [0,1]:
          blockMeshDict.write("    hex ("+" ".join([str(int(4*k+j+x+12*i)%(12*z)+16) for x in [0,1,13,12,4,5,17,16]])+") rotor"+(str(i) if (k is 0 and j is 1) else str())+" ("+str(Noxo)+" "+str(Not)+" "+str(Nori if (k is 0) else Noro)+") simpleGrading (1 1 1)\n");
    blockMeshDict.write(");\n");
    blockMeshDict.write("\n");
    blockMeshDict.write("edges\n");
    blockMeshDict.write("(\n");
    for i in range(4):
        alpha=0.5*np.pi*(0.5+i)
        blockMeshDict.write("    arc "+str((i*4)%16)+" "+str((4*(i+1))%16)+" ("+str(-0.5*L)+" "+str(Rs*np.sin(alpha))+" "+str(Rs*np.cos(alpha))+")\n");
        blockMeshDict.write("    arc "+str((i*4+1)%16)+" "+str((4*(i+1)+1)%16)+" ("+str(0.5*L)+" "+str(Rs*np.sin(alpha))+" "+str(Rs*np.cos(alpha))+")\n");
        blockMeshDict.write("    arc "+str((i*4+3)%16)+" "+str((4*(i+1)+3)%16)+" ("+str(-0.5*L)+" "+str(Ro*np.sin(alpha))+" "+str(Ro*np.cos(alpha))+")\n");
        blockMeshDict.write("    arc "+str((i*4+2)%16)+" "+str((4*(i+1)+2)%16)+" ("+str(0.5*L)+" "+str(Ro*np.sin(alpha))+" "+str(Ro*np.cos(alpha))+")\n");
    for i in range(z):
      alpha=2.0*np.pi/z*(0.5+i)
      for j, r in enumerate([Rs, Rb,Rb+Do]):
        r
        for k, x in enumerate([-Do-0.5*Lb, -0.5*Lb, 0.5*Lb, Do+0.5*Lb]):
          xcoord=str(x)
          blockMeshDict.write("    arc "+str(int(12*i+4*j+k)%(12*z)+16)+" "+str(int(12*(i+1)+4*j+k)%(12*z)+16)+" ("+str(xcoord)+" "+str(r*np.sin(alpha))+" "+str(r*np.cos(alpha))+")\n");
    blockMeshDict.write(");\n");
    blockMeshDict.write("\n");
    blockMeshDict.write("boundary\n");
    blockMeshDict.write("(\n");
    blockMeshDict.write("    inlet\n");
    blockMeshDict.write("    {\n");
    blockMeshDict.write("        type patch;\n");
    blockMeshDict.write("        faces\n");
    blockMeshDict.write("        (\n");
    for i in range(4):
        blockMeshDict.write("            ("+" ".join([str(int(x+4*i)%16) for x in [0,3,7,4]])+")\n");
    blockMeshDict.write("        );\n");
    blockMeshDict.write("    }\n");
    blockMeshDict.write("    outlet\n");
    blockMeshDict.write("    {\n");
    blockMeshDict.write("        type patch;\n");
    blockMeshDict.write("        faces\n");
    blockMeshDict.write("        (\n");
    for i in range(4):
        blockMeshDict.write("            ("+" ".join([str(int(x+4*i)%16) for x in [5,6,2,1]])+")\n");
    blockMeshDict.write("        );\n");
    blockMeshDict.write("    }\n");
    blockMeshDict.write("    fixedWall\n");
    blockMeshDict.write("    {\n");
    blockMeshDict.write("        type wall;\n");
    blockMeshDict.write("        faces\n");
    blockMeshDict.write("        (\n");
    for i in range(4):
        blockMeshDict.write("            ("+" ".join([str(int(x+4*i)%16) for x in [0,4,5,1]])+")\n");
        blockMeshDict.write("            ("+" ".join([str(int(x+4*i)%16) for x in [3,2,6,7]])+")\n");
    blockMeshDict.write("        );\n");
    blockMeshDict.write("    }\n");
    blockMeshDict.write("    nonRotatingWall\n");
    blockMeshDict.write("    {\n");
    blockMeshDict.write("        type wall;\n");
    blockMeshDict.write("        faces\n");
    blockMeshDict.write("        (\n");
    for i in range(z):
      for j in range(2):
        blockMeshDict.write("            ("+" ".join([str(int(12*i+2*j+x)%(12*z)+16) for x in [12,13,1,0]])+")\n");
    blockMeshDict.write("        );\n");
    blockMeshDict.write("    }\n");
    blockMeshDict.write("    rotatingWall\n");
    blockMeshDict.write("    {\n");
    blockMeshDict.write("        type wall;\n");
    blockMeshDict.write("        faces\n");
    blockMeshDict.write("        (\n");
    for i in range(z):
      blockMeshDict.write("            ("+" ".join([str(int(12*i+x)%(12*z)+16) for x in [13,14,2,1]])+")\n");
    blockMeshDict.write("        );\n");
    blockMeshDict.write("    }\n");
    blockMeshDict.write("    inner\n");
    blockMeshDict.write("    {\n");
    blockMeshDict.write("        type bellerophon;\n");
    blockMeshDict.write("        donorZone stator;\n");
    blockMeshDict.write("        faces\n");
    blockMeshDict.write("        (\n");
    for i in range(z):
      for j in range(2):
        blockMeshDict.write("            ("+" ".join([str(int(12*i+j*4+x)%(12*z)+16) for x in [0,4,16,12]])+")\n");
        blockMeshDict.write("            ("+" ".join([str(int(12*i+j*4+x)%(12*z)+16) for x in [15,19,7,3]])+")\n");
      for j in range(3):
        blockMeshDict.write("            ("+" ".join([str(int(12*i+j+x)%(12*z)+16) for x in [8,9,21,20]])+")\n");
    blockMeshDict.write("        );\n");
    blockMeshDict.write("    }\n");
    blockMeshDict.write(");\n");
    blockMeshDict.write("\n");
    blockMeshDict.write("// ************************************************************************* //\n");
    blockMeshDict.write("\n");
    blockMeshDict.close();
    
    os.system("blockMesh")
    
    topoSetDict = open("system/topoSetDict","w")
    topoSetDict.write("/*--------------------------------*- C++ -*----------------------------------*\\\n");
    topoSetDict.write("| =========                 |                                                 |\n");
    topoSetDict.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n");
    topoSetDict.write("|  \\\\    /   O peration     | Version:  2.3.0                                 |\n");
    topoSetDict.write("|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n");
    topoSetDict.write("|    \\\\/     M anipulation  |                                                 |\n");
    topoSetDict.write("\\*---------------------------------------------------------------------------*/\n");
    topoSetDict.write("FoamFile\n");
    topoSetDict.write("{\n");
    topoSetDict.write("    version     2.0;\n");
    topoSetDict.write("    format      ascii;\n");
    topoSetDict.write("    class       dictionary;\n");
    topoSetDict.write("    object      topoSetDict;\n");
    topoSetDict.write("}\n");
    topoSetDict.write("\n");
    topoSetDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n");
    topoSetDict.write("\n");
    topoSetDict.write("actions\n");
    topoSetDict.write("(\n");
    for i in range(z):
        topoSetDict.write("    {\n");
        topoSetDict.write("        name    cs_"+str(i)+";\n");
        topoSetDict.write("        type    cellSet;\n");
        topoSetDict.write("        action  new;\n");
        topoSetDict.write("        source  zoneToCell;\n");
        topoSetDict.write("        sourceInfo\n");
        topoSetDict.write("        {\n");
        topoSetDict.write("            name rotor"+str(i)+";\n");
        topoSetDict.write("        }\n");
        topoSetDict.write("    }\n");
        topoSetDict.write("\n");
    for i in range(z):
        j=(i+1)%z
        topoSetDict.write("    {\n");
        topoSetDict.write("        name    set_"+str(i)+"_"+str(j)+";\n");
        topoSetDict.write("        type    faceSet;\n");
        topoSetDict.write("        action  new;\n");
        topoSetDict.write("        source  cellToFace;\n");
        topoSetDict.write("        sourceInfo\n");
        topoSetDict.write("        {\n");
        topoSetDict.write("            set cs_"+str(i)+";\n");
        topoSetDict.write("            option all;\n");
        topoSetDict.write("        }\n");
        topoSetDict.write("    }\n");
        topoSetDict.write("\n");
        topoSetDict.write("    {\n");
        topoSetDict.write("        name    set_"+str(i)+"_"+str(j)+";\n");
        topoSetDict.write("        type    faceSet;\n");
        topoSetDict.write("        action  subset;\n");
        topoSetDict.write("        source  cellToFace;\n");
        topoSetDict.write("        sourceInfo\n");
        topoSetDict.write("        {\n");
        topoSetDict.write("            set cs_"+str(j)+";\n");
        topoSetDict.write("            option all;\n");
        topoSetDict.write("        }\n");
        topoSetDict.write("    }\n");
        topoSetDict.write("\n");
    for i in range(z):
        j=(i+1)%z
        topoSetDict.write("    {\n");
        topoSetDict.write("        name    buffleZone;\n");
        topoSetDict.write("        type    faceZoneSet;\n");
        if i is 0:
            topoSetDict.write("        action  new;\n");
        else:
            topoSetDict.write("        action  add;\n");
        topoSetDict.write("        source  setToFaceZone;\n");
        topoSetDict.write("        sourceInfo\n");
        topoSetDict.write("        {\n");
        topoSetDict.write("            faceSet set_"+str(i)+"_"+str(j)+";\n");
        topoSetDict.write("        }\n");
        topoSetDict.write("    }\n");
        topoSetDict.write("\n");
    for i in range(z):
        topoSetDict.write("    {\n");
        topoSetDict.write("        name    rotor"+str(i)+";\n");
        topoSetDict.write("        type    cellZoneSet;\n");
        topoSetDict.write("        action  remove;\n");
        topoSetDict.write("    }\n");
        topoSetDict.write("\n");
    topoSetDict.write("    {\n");
    topoSetDict.write("        name    rotor;\n");
    topoSetDict.write("        type    cellZoneSet;\n");
    topoSetDict.write("        action  clear;\n");
    topoSetDict.write("    }\n");
    topoSetDict.write("\n");
    topoSetDict.write("    {\n");
    topoSetDict.write("        name    rotorSet;\n");
    topoSetDict.write("        type    cellSet;\n");
    topoSetDict.write("        action  new;\n");
    topoSetDict.write("        source  zoneToCell;\n");
    topoSetDict.write("        sourceInfo\n");
    topoSetDict.write("        {\n");
    topoSetDict.write("            name stator;\n");
    topoSetDict.write("        }\n");
    topoSetDict.write("    }\n");
    topoSetDict.write("\n");
    topoSetDict.write("    {\n");
    topoSetDict.write("        name    rotorSet;\n");
    topoSetDict.write("        type    cellSet;\n");
    topoSetDict.write("        action  invert;\n");
    topoSetDict.write("    }\n");
    topoSetDict.write("\n");
    topoSetDict.write("    {\n");
    topoSetDict.write("        name    rotor;\n");
    topoSetDict.write("        type    cellZoneSet;\n");
    topoSetDict.write("        action  new;\n");
    topoSetDict.write("        source  setToCellZone;\n");
    topoSetDict.write("        sourceInfo\n");
    topoSetDict.write("        {\n");
    topoSetDict.write("            set rotorSet;\n");
    topoSetDict.write("        }\n");
    topoSetDict.write("    }\n");
    topoSetDict.write("\n");
    topoSetDict.write(");\n");
    topoSetDict.write("\n");
    topoSetDict.write("// ************************************************************************* //\n");
    topoSetDict.write("\n");
    topoSetDict.close();
    
    os.system("topoSet")
    
    createBafflesDict = open("system/createBafflesDict","w")
    createBafflesDict.write("/*--------------------------------*- C++ -*----------------------------------*\\\n");
    createBafflesDict.write("| =========                 |                                                 |\n");
    createBafflesDict.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n");
    createBafflesDict.write("|  \\\\    /   O peration     | Version:  dev                                   |\n");
    createBafflesDict.write("|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n");
    createBafflesDict.write("|    \\\\/     M anipulation  |                                                 |\n");
    createBafflesDict.write("\\*---------------------------------------------------------------------------*/\n");
    createBafflesDict.write("FoamFile\n");
    createBafflesDict.write("{\n");
    createBafflesDict.write("    version     2.0;\n");
    createBafflesDict.write("    format      ascii;\n");
    createBafflesDict.write("    class       dictionary;\n");
    createBafflesDict.write("    object      createBafflesDict;\n");
    createBafflesDict.write("}\n");
    createBafflesDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n");
    createBafflesDict.write("\n");
    createBafflesDict.write("// Whether to convert internal faces only (so leave boundary faces intact).\n");
    createBafflesDict.write("// This is only relevant if your face selection type can pick up boundary\n");
    createBafflesDict.write("// faces.\n");
    createBafflesDict.write("internalFacesOnly true;\n");
    createBafflesDict.write("\n");
    createBafflesDict.write("\n");
    createBafflesDict.write("// Baffles to create.\n");
    createBafflesDict.write("baffles\n");
    createBafflesDict.write("{\n");
    createBafflesDict.write("    baffleFaces\n");
    createBafflesDict.write("    {\n");
    createBafflesDict.write("        //- Use predefined faceZone to select faces and orientation.\n");
    createBafflesDict.write("        type        faceZone;\n");
    createBafflesDict.write("        zoneName    buffleZone;\n");
    createBafflesDict.write("\n");
    createBafflesDict.write("        patchPairs\n");
    createBafflesDict.write("        {\n");
    createBafflesDict.write("            type            wall;\n");
    createBafflesDict.write("\n");
    createBafflesDict.write("            patchFields\n");
    createBafflesDict.write("            {\n");
    createBafflesDict.write("            }\n");
    createBafflesDict.write("        }\n");
    createBafflesDict.write("    }\n");
    createBafflesDict.write("}\n");
    createBafflesDict.write("\n");
    createBafflesDict.write("\n");
    createBafflesDict.write("// ************************************************************************* //\n");
    createBafflesDict.close()
    
    os.system("createBaffles -overwrite")
    
    os.system("rm -rf constant/polyMesh/sets")
    
    topoSetDict = open("system/topoSetDict","w")
    topoSetDict.write("/*--------------------------------*- C++ -*----------------------------------*\\\n");
    topoSetDict.write("| =========                 |                                                 |\n");
    topoSetDict.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n");
    topoSetDict.write("|  \\\\    /   O peration     | Version:  2.3.0                                 |\n");
    topoSetDict.write("|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n");
    topoSetDict.write("|    \\\\/     M anipulation  |                                                 |\n");
    topoSetDict.write("\\*---------------------------------------------------------------------------*/\n");
    topoSetDict.write("FoamFile\n");
    topoSetDict.write("{\n");
    topoSetDict.write("    version     2.0;\n");
    topoSetDict.write("    format      ascii;\n");
    topoSetDict.write("    class       dictionary;\n");
    topoSetDict.write("    object      topoSetDict;\n");
    topoSetDict.write("}\n");
    topoSetDict.write("\n");
    topoSetDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n");
    topoSetDict.write("\n");
    topoSetDict.write("actions\n");
    topoSetDict.write("(\n");
    topoSetDict.write("    {\n");
    topoSetDict.write("        name    holeSource;\n");
    topoSetDict.write("        type    cellSet;\n");
    topoSetDict.write("        action  new;\n");
    topoSetDict.write("        source  cylinderToCell;\n");
    topoSetDict.write("        sourceInfo\n");
    topoSetDict.write("        {\n");
    topoSetDict.write("            p1 ("+str(-Do-0.5*Lb)+" 0 0);\n");
    topoSetDict.write("            p2 ("+str(Do+0.5*Lb)+" 0 0);\n");
    topoSetDict.write("            radius "+str(Rb+Do)+";\n");
    topoSetDict.write("        }\n");
    topoSetDict.write("    }\n");
    topoSetDict.write("\n");
    topoSetDict.write("    {\n");
    topoSetDict.write("        name    holeSource;\n");
    topoSetDict.write("        type    cellSet;\n");
    topoSetDict.write("        action  delete;\n");
    topoSetDict.write("        source  zoneToCell;\n");
    topoSetDict.write("        sourceInfo\n");
    topoSetDict.write("        {\n");
    topoSetDict.write("            name rotor;\n");
    topoSetDict.write("        }\n");
    topoSetDict.write("    }\n");
    topoSetDict.write(");\n");
    topoSetDict.write("\n");
    topoSetDict.write("// ************************************************************************* //\n");
    topoSetDict.write("\n");
    topoSetDict.close();
    
    os.system("topoSet")
    os.system("holeCutter")
    os.system("subsetMesh liveCells -overwrite")
    
    createPatchDict = open("system/createPatchDict", "w")
    createPatchDict.write("/*--------------------------------*- C++ -*----------------------------------*\\\n");
    createPatchDict.write("| =========                 |                                                 |\n");
    createPatchDict.write("| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n");
    createPatchDict.write("|  \\\\    /   O peration     | Version:  2.3.x                                 |\n");
    createPatchDict.write("|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |\n");
    createPatchDict.write("|    \\\\/     M anipulation  |                                                 |\n");
    createPatchDict.write("\\*---------------------------------------------------------------------------*/\n");
    createPatchDict.write("FoamFile\n");
    createPatchDict.write("{\n");
    createPatchDict.write("    version     2.0;\n");
    createPatchDict.write("    format      ascii;\n");
    createPatchDict.write("    class       dictionary;\n");
    createPatchDict.write("    object      createPatchDict;\n");
    createPatchDict.write("}\n");
    createPatchDict.write("\n");
    createPatchDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n");
    createPatchDict.write("pointSync false;\n");
    createPatchDict.write("\n");
    createPatchDict.write("patches\n");
    createPatchDict.write("(\n");
    createPatchDict.write("{\n");
    createPatchDict.write("    name outer;\n");
    createPatchDict.write("    patchInfo\n");
    createPatchDict.write("    {\n");
    createPatchDict.write("        type bellerophon;\n");
    createPatchDict.write("        donorZone rotor;\n");
    createPatchDict.write("    }\n");
    createPatchDict.write("    constructFrom patches;\n");
    createPatchDict.write("    patches (oldInternalFaces);\n");
    createPatchDict.write("}\n");
    createPatchDict.write("{\n");
    createPatchDict.write("    name rotatingWall;\n");
    createPatchDict.write("    patchInfo\n");
    createPatchDict.write("    {\n");
    createPatchDict.write("        type wall;\n");
    createPatchDict.write("    }\n");
    createPatchDict.write("    constructFrom patches;\n");
    createPatchDict.write("    patches (rotatingWall baffleFaces_slave baffleFaces_master);\n");
    createPatchDict.write("}\n");
    createPatchDict.write(");\n");
    createPatchDict.write("\n");
    createPatchDict.write("// ************************************************************************* //\n");
    createPatchDict.close()
    
    os.system("createPatch -overwrite")

    os.system("rm -rf 0")
    os.system("cp -r 0.org 0")

if __name__ == "__main__":
    import sys
    z=False
    if len(sys.argv) > 1:
        try:
            z=int(sys.argv[1])
        except:
            pass
    print(z)
    if z is not False:
        makeMesh(z)
    else:
        makeMesh()
