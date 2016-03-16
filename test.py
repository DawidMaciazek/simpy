from simpy import analyze

tr = analyze.traj("test/test.lammpstrj")
tr.read(1)

