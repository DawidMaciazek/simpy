from simpy import analyze
from simpy import write
import numpy as np

traj = analyze.Traj("se_01107.lammpstrj")

frame = traj.read()

print frame.keys()

frame.remove_element('Kr')
frame.remove_lesser_than(-30,2)

print frame['box']
frame.periodic([ [-1,1], [-1,0], [0,0]])
print frame['box']
