from simpy import analyze
from simpy import write
from os import listdir
import numpy as np
import sys

def calc_displacment(dir):
    traj = analyze.Traj(dir)

    # remove new ion form start and from end if implanted
    sframe = traj.read()
    scoord = sframe['coord'][:-1]
    projectile_id = sframe['id'][-1]
    sid = sframe['id'][:-1]

    eframe = traj.read()
    if eframe['id'][-1] == projectile_id:
        ecoord = eframe['coord'][:-1]
        eid = eframe['id'][:-1]
        etype = eframe['type'][:-1]
    else:
        ecoord = eframe['coord']
        eid = eframe['id']
        etype = eframe['type']


    shift = np.empty(ecoord.shape, dtype=float)
    shift_si = np.empty(ecoord.shape, dtype=bool)

    iend = 0
    for i in range(len(sid)):
        if sid[i] != eid[iend]:
            continue

        shift[iend] = ecoord[iend] - scoord[i]
        shift_si[iend] = (etype[iend] == 2)
        iend += 1

    if iend != len(ecoord):
        print "ERROR !!!!"
        print "missing atoms"
        sys.exit()
    #shift = shift[:,0]

    print sum(shift[:,0])
    print sum(shift[:,1])
    print sum(shift[:,2])

    # RMS and surf lvl analyze
    # ------------------------
    sframe.remove_type(1)
    sframe.remove_lesser_than(-30, 2)
    old_box = sframe['box']

    sframe.periodic([[-1,1], [-1,1], [0,0]])
    rms = analyze.RMS(sframe['coord'], 2.1, 2.1, old_box)
    print rms.compute()

    name = dir.split('/')[-1].split(".")[0] + ".xyz"
    write.xyz(name).write(rms.get_frame())



def get_rms():
    frame.remove_type(1)
    frame.remove_lesser_than(-30,2)
    old_box = frame['box']

    frame.periodic([[-1,1], [-1,1], [0,0]])
    rms = analyze.RMS(frame['coord'], 2.1, 2.1, old_box)

    rms.compute()
    frame = rms.get_frame()
    write.xyz("out.xyz").write(frame)



basedir = "/home/dawid/dumps/hobler_redist/si_kr_2kev/dumpse"
dirs = sorted(listdir(basedir))


stop = 0
for dir in dirs:

    dirr = basedir + "/" + dir
    calc_displacment(dirr)
    if stop == 10:
        break
    stop += 1

