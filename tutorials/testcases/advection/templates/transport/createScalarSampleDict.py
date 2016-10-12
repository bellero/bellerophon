#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 12:44:29 2015

@author: Andreas Groß

Create sampleDict for sampling T in plane
define path as argument is optional
"""

import math
import sys
import os

def nFromBoundAndDelta(r,d):
  """
  Calculate number of steps from range and delta
  """
  return int(math.ceil((r[1]-r[0])/d)+1.0)

def xFromRangeNAndI(r,n,i):
  """
  Calculate current coordinate from range, spacing and iterator
  """
  if int(n) == 1:
    return r[0]
  else:
    return r[0]+(r[1]-r[0])/float(n-1)*float(i)

delta = 0.005

x_range = [0.0, 2.5]                  # Domain boundary in x 
y_range = [0.0, 1.0]                  # Domain boundary in y
z_range = [0.05, 0.05]                # Domain boundary in z

offset = [0.0002,0.0002,0]           # offset from calculated points
                                      # shifted, if resulting point outside of domain
         
path = None
if len(sys.argv) > 1:
  path = sys.argv[1]

nx=nFromBoundAndDelta(x_range,delta)
ny=nFromBoundAndDelta(y_range,delta)
nz=nFromBoundAndDelta(z_range,delta)

print("Size: "+str(nx)+" x "+str(ny)+" x "+str(nz))

fileName = "sampleDict"
if path:
  if os.path.isdir(path):
    fileName = path+"/"+fileName
  else:
    print("specified path not found!")
    raise IOError
    

sampleDict = open(fileName,'w')

sampleDict.write("/*--------------------------------*- C++ -*----------------------------------*\\\n")
sampleDict.write("| =========                 |                                                 |\n")
sampleDict.write("| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n")
sampleDict.write("|  \\    /   O peration     | Version:  2.3.0                                 |\n")
sampleDict.write("|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |\n")
sampleDict.write("|    \\/     M anipulation  |                                                 |\n")
sampleDict.write("\*---------------------------------------------------------------------------*/\n")
sampleDict.write("FoamFile\n")
sampleDict.write("{\n")
sampleDict.write("    version     2.0;\n")
sampleDict.write("    format      ascii;\n")
sampleDict.write("    class       dictionary;\n")
sampleDict.write("    location    \"system\";\n")
sampleDict.write("    object      sampleDict;\n")
sampleDict.write("}\n")
sampleDict.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n")
sampleDict.write("\n")
sampleDict.write("interpolationScheme cellPoint;\n")
sampleDict.write("\n")
sampleDict.write("setFormat       raw;\n")
sampleDict.write("\n")
sampleDict.write("sets\n")
sampleDict.write("(\n")
sampleDict.write("    sampleCloud\n")
sampleDict.write("    {\n")
sampleDict.write("        type    cloud;\n")
sampleDict.write("        axis    xyz;\n")
sampleDict.write("        points\n")
sampleDict.write("        (\n")

for i in xrange(nx):
  x = xFromRangeNAndI(x_range,nx,i)
  for j in xrange(ny):
    y = xFromRangeNAndI(y_range,ny,j)
    for k in xrange(nz):
      z = xFromRangeNAndI(z_range,nz,k)
      if (x+offset[0])>x_range[1]:
        xcorr=x-offset[0]
      else:
        xcorr=x+offset[0]
      if (y+offset[1])>y_range[1]:
        ycorr=y-offset[1]
      else:
        ycorr=y+offset[1]
      if (z+offset[2])>z_range[1]:
        zcorr=z-offset[2]
      else:
        zcorr=z+offset[2]
      
      sampleDict.write("            ("+str(xcorr)+" "+str(ycorr)+" "+str(zcorr)+")\n")

sampleDict.write("        );\n")
sampleDict.write("    }\n")
sampleDict.write(");\n")
sampleDict.write("\n")
sampleDict.write("fields          ( T );\n")
sampleDict.write("\n")
sampleDict.write("\n")
sampleDict.write("// ************************************************************************* //\n")
sampleDict.close()