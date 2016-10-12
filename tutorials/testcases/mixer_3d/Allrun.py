#!/usr/bin/env/python

from makeMesh import makeMesh

makeMesh()

import os
import shutil

if os.path.isdir("0"):
    shutil.rmtree("0")
shutil.copytree("0.org","0")


os.system("decomposePar")
os.system("mpirun -np 4 simpleFoam -parallel")
os.system("reconstructPar -latestTime")

