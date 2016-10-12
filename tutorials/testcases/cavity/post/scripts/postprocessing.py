#!/usr/bin/python3

dimension = 0.1
nCells = 30
dX = dimension/nCells

z = 0.005

f = open('../results/probeLocations', 'w')

for i in range(0, nCells):
  x = (i+0.5)*dX
  for j in range(0, nCells):
    y = (j+0.5)*dX
    f.write('    ({0} {1} {2})\n'.format(x, y, z))

f.close()
