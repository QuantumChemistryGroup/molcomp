#!/usr/bin/python
import sys, random

f = open(sys.argv[1], 'r')
coord = f.readlines() # read every line
f.close()
d = str.split(coord[0])
if len(d) == 1:
  MOL = True
else: MOL = False 

if MOL:
  coord = coord[2:]

coordinates = []
for i in coord:
  i = str.split(i)
  if len(i) > 3:
    coordinates.append(i)

for i in range(len(coordinates)):
  coordinates[i][3] = str(float(coordinates[i][3])*-1.0)
random.seed(42)
random.shuffle(coordinates)

#printing
f = open(sys.argv[1][0:-4]+'mirror.xyz', 'w')
f.write(str(len(coordinates))+'\n')
f.write('Molecule\n')
for i in coordinates:
  f.write(i[0]+'   '+i[1]+'   '+i[2]+ '   '+i[3]+'\n')
f.close()
