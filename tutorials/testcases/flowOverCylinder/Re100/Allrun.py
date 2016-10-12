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

import caseSetup
import os
import shutil
#from multiprocessing import Pool

if __name__ == '__main__':
  s = caseSetup.Setup()

  if not os.path.exists(s.results):
    os.makedirs(s.results)
  
  if os.path.exists(s.calculations):
    shutil.rmtree(s.calculations)
  
  os.mkdir(s.calculations)

  caseSetup.doSingle(s)
  caseSetup.doOverset(s)

  print("Done with calculations, post...")
  caseSetup.doPostprocessing(s)
