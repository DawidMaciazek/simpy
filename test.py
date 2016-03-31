from simpy import analyze
from simpy import write

tr = analyze.traj("test/surf_test.lammpstrj")
frames = tr.read(-1)


print frames[0][0]["format"]


w = write.xyz("out.xyz")

w.write(frames[0][2], frames[0][3])
w.write(frames[1][2], frames[1][3])

w.close()
