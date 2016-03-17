from simpy import analyze

tr = analyze.traj("test/test.lammpstrj")
frame = tr.read(-1)

