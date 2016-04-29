from simpy import analyze
from simpy import write


def single(types, coords):
    wr = write.xyz('single.xyz')
    wr.write(types, coords)
    wr.close()



traj = analyze.Traj('../resources/surf_test.xyz', 'xyz')
frames = traj.read(1)

# test single frame
single(frames[0][1], frames[0][2])
