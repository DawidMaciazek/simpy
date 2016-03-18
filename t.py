import numpy

f = open('out.xyz', 'r')

atoms = int(f.readline())
f.readline()

z_array = numpy.zeros(atoms)
for i in xrange(atoms):
    l = f.readline()
    s = l.split()
    z_array[i] = float(s[3])

print z_array[0:20]
ave =numpy.average(z_array)
print ave

summ=0
for i in xrange(atoms):
    res=ave-z_array[i]
    summ+=res*res

import math
print math.sqrt(summ/atoms)
