#!/usr/bin/env python
  
import caseSetup
import shutil
import os
from multiprocessing import Pool

if __name__ == '__main__':
  s = caseSetup.Setup()

  doPar = False

  if not os.path.exists(s.postDir):
    os.mkdir(s.postDir)
  
  #if os.path.exists(s.resultDir):
  #  shutil.rmtree(s.resultDir)
  #os.mkdir(s.resultDir)
  
  #clean all old calculations
  #if(os.path.exists(s.calculationsDir)):
  #  shutil.rmtree(s.calculationsDir)
  
  if doPar:
    p = Pool(maxtasksperchild=1)

    job=[]
    job.append(p.map_async(caseSetup.doSingle,[s]))
    job.append(p.map_async(caseSetup.doPerfectMatching,[s]))
    job.append(p.map_async(caseSetup.doShifted,[s]))
    job.append(p.map_async(caseSetup.doTri,[s]))
    job.append(p.map_async(caseSetup.doPoly,[s]))
    p.close()
    p.join()

  else:
    #caseSetup.doSingle(s)
    #caseSetup.doPerfectMatching(s)
    #caseSetup.doShifted(s)
    #caseSetup.doTri(s)
    caseSetup.doPoly(s)
