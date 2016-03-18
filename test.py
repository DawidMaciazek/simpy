from simpy import analyze

tr = analyze.traj("test/surf_test.lammpstrj")
frames = tr.read(-1)

frame = frames[1]


rms = analyze.rms(frame[3], 2, 1.5)
import numpy

x = numpy.arange(-57.5, 58, 3)
y = numpy.arange(-57.5, 58, 3)
X, Y =numpy.meshgrid(x, y)

probe_array = numpy.array([X.flatten(), Y.flatten()]).T

out_array = rms.get_surf_array_oa2d(probe_array)

print len(out_array)


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


