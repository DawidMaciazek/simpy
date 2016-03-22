from simpy import analyze

tr = analyze.traj("test/surf_test.lammpstrj")
frames = tr.read(-1)

frame = frames[1]

import numpy
for i in xrange(len(frame[2])):
    if frame[2][i] == 'Ga':
        #print "galll"
        frame[3][i][2] = -20.00


rms = analyze.rms(frame[3], 1, 1)

x = numpy.arange(-57.5, 58, 2.9)
y = numpy.arange(-57.5, 58, 2.9)
X, Y =numpy.meshgrid(x, y)

probe_array = numpy.array([X.flatten(), Y.flatten()]).T

out_array = rms.get_surf_array_oa2d(probe_array, lc_rep=3)


from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
xxx = out_array[:,0]
yyy = out_array[:,1]
zzz = out_array[:,2]
print xxx[:5]
print yyy[:5]
print zzz[:5]
print out_array[:5]

ax.scatter(xxx, yyy, zzz, c='r')

ax.set_xlabel("X label")
ax.set_xlabel("Y label")
ax.set_xlabel("Z label")

print ax.set_zlim([-50, 50])

plt.show()



#ff_out = numpy.fft.rfft2(out_array)







print len(out_array)

print " OUTTT"
ave = numpy.average(out_array[:,2])

summ=0
for i in xrange(len(out_array)):
    res = ave-out_array[i][2]
    summ += res*res

print summ

import math
print "RMS:", math.sqrt(summ/len(out_array))

#lc = rms.logical_cells











n = len(out_array)


for i in xrange(n):
    a = out_array[i]
    if numpy.isnan(a[2]):
        print "NONE"
        out_array[i][2]=-100

f = open("out.xyz", 'w')
f.write("%i\n%i" % (n, n))




for i in xrange(n):
    a = out_array[i]
    f.write("\nH %f %f %f" % (a[0], a[1], a[2]))

