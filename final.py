from simpy import analyze
import numpy
import math

tr = analyze.traj("/home/dawid/dumps/silicon/lj_full_100.lammpstrj")

x = numpy.arange(-57.5, 58, 1.0)
y = numpy.arange(-57.5, 58, 1.0)
X, Y =numpy.meshgrid(x, y)
probe_array = numpy.array([X.flatten(), Y.flatten()]).T


frame = tr.read(1)

output = open("out_100_rms.txt", 'w')


frame_cnt = 0
while frame:
    frame = frame[0]
    frame_cnt += 1
    for i in xrange(len(frame[2])):
        if frame[2][i] == 'Ga':
            frame[3][2] = -10.0

    rms = analyze.rms(frame[3], 1, 1)
    out_array = rms.get_surf_array_oa2d(probe_array, lc_rep=3)

    ave = numpy.average(out_array[:,2])
    summ = 0
    for i in xrange(len(out_array)):
        res = ave-out_array[i][2]
        summ += res*res

    print "%i %f" % (frame_cnt, math.sqrt(summ/len(out_array)))
    output.write("%i %f\n" % (frame_cnt, math.sqrt(summ/len(out_array))))


    frame = tr.read(1)

output.close()
