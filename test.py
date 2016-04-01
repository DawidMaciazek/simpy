from simpy import analyze
from simpy import write

tr = analyze.traj("test/surf_test.lammpstrj")
frames = tr.read(-1)


print frames[0][0]["format"]

print frames[0][1].shape

w = write.xyz("out.xyz")
new_types = frames[0][2]
w.write(new_types, frames[0][3], new_types)

w.close()
