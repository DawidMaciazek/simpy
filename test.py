from simpy import analyze

tr = analyze.traj("test/surf_test.lammpstrj")
frames = tr.read(-1)

frame = frames[1]


rms = analyze.rms(frame[3], 1, 1)

lc = rms.logical_cells

f = open("out.xyz", 'w')
n = len(lc)*len(lc[0])
f.write("%i\n%i" % (n, n))
for i in xrange(len(lc)):
    for j in xrange(len(lc[0])):
        f.write("\nH %f %f %f" % (i, j, lc[i][j]))


