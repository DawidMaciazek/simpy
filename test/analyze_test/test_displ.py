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
    stype = sframe['type'][:-1]


    eframe = traj.read()
    projcoord = sframe['coord'][-1]
    if eframe['id'][-1] == projectile_id:
        ecoord = eframe['coord'][:-1]
        eid = eframe['id'][:-1]
        etype = eframe['type'][:-1]
    else:
        ecoord = eframe['coord']
        eid = eframe['id']
        etype = eframe['type']

    start_z = np.empty(len(ecoord), dtype=float)


    shift = np.empty(ecoord.shape, dtype=float)
    shift_si_flag = np.empty(ecoord.shape[0], dtype=bool)

    iend = 0
    sputtered = []
    sputtered_si_flag = []
    for i in range(len(sid)):
        if len(eid) == iend:
            sputtered.append(scoord[i])
            sputtered_si_flag.append(stype[i] == 2)
            continue

        if sid[i] != eid[iend]:
            sputtered.append(scoord[i])
            sputtered_si_flag.append(stype[i] == 2)
            continue

        shift[iend] = ecoord[iend] - scoord[i]
        shift_si_flag[iend] = (etype[iend] == 2)

        start_z[iend] = scoord[i][2]

        iend += 1

    sputtered = np.array(sputtered)
    sputtered_si_flag = np.array(sputtered_si_flag)

    # sputtered atoms ....
    if len(sputtered_si_flag) != 0:
        #  (Si sputt, Kr sputt)
        sputt_list = (sputtered[sputtered_si_flag], sputtered[np.logical_not(sputtered_si_flag)])
    else:

        sputt_list = (None, None)

    # remove atoms below cutoff
    cutoff = 0.01
    above_cutoff = np.abs(shift[:,2]) >= cutoff

    shifted_list = (shift[np.logical_and(above_cutoff, shift_si_flag)],
              shift[np.logical_and(above_cutoff, np.logical_not(shift_si_flag))])

    start_z_list = (start_z[np.logical_and(above_cutoff, shift_si_flag)],
              start_z[np.logical_and(above_cutoff, np.logical_not(shift_si_flag))])
    print len(shifted_list[0]), len(shifted_list[1]), " <<<<<<<< "

    # projectile
    # si atoms


    if iend != len(ecoord):
        print "ERROR !!!!"
        print "missing atoms"
        print dir
        sys.exit()

    print sputtered
    print sum(shift[:,0])
    print sum(shift[:,1])
    print sum(shift[:,2])


    # ------------------------
    # RMS and surf lvl analyze
    # ------------------------
    sframe.remove_type(1)
    sframe.remove_lesser_than(-30, 2)
    old_box = sframe['box']

    sframe.periodic([[-1,1], [-1,1], [0,0]])
    rms = analyze.RMS(sframe['coord'], 2.1, 2.1, old_box)

    ave_surf, rms_val = rms.compute()

    return (sputt_list, shifted_list, start_z_list, ave_surf, rms_val, projcoord)


    bin_size = 2
    for i in range(-18,5)[::-1]:
        start = ave_surf + i*bin_size
        stop = start - bin_size

        layer = np.logical_and(start_z <= start, start_z > stop)
        print "l ({}, {}) : {}   // atoms = {}".format(start, stop, sum(shift[layer][:,2]), sum(layer))






    #name = dir.split('/')[-1].split(".")[0] + ".xyz"
    #write.xyz(name).write(rms.get_frame())



def get_rms():
    frame.remove_type(1)
    frame.remove_lesser_than(-30,2)
    old_box = frame['box']

    frame.periodic([[-1,1], [-1,1], [0,0]])
    rms = analyze.RMS(frame['coord'], 2.1, 2.1, old_box)

    rms.compute()
    frame = rms.get_frame()
    write.xyz("out.xyz").write(frame)


def write_redist(traj_num, redist, start_z_list, outfile, ave_surf, rms_val):
    # redistributive
    outfile.write("# {} {} {} {}\n".format(traj_num, len(redist), ave_surf, rms_val))

    for i in xrange(len(redist)):
        outfile.write("{}  {} {} {}\n".format(start_z_list[i] - ave_surf, redist[i][0], redist[i][1],redist[i][2]))


def write_sputt(traj_num, sputt, outfile, ave_surf, rms_val, projcoord):
    ## box size ...

    box_half_len = 54.31
    box_len = box_half_len*2

    backward_f = 0.25
    backward_a = backward_f * box_len

    z_proj_path = projcoord[2] - ave_surf
    x_proj_path = np.sqrt(3) * z_proj_path

    x_hit_position = (x_proj_path + projcoord[0] + box_half_len) % box_len
    x_shift = backward_a - x_hit_position

    y_hit_position = (projcoord[1] + box_half_len) % box_len
    y_shift = box_half_len - y_hit_position

    print  "---------------------"
    print "x_shift:", x_shift
    print "y_shift:", y_shift
    print ave_surf
    print projcoord
    print x_hit_position
    print  "---------------------"

    if sputt is None:
        print "skipping sputtering - None"
        outfile.write("# {} {}\n".format(traj_num, 0, ave_surf, rms_val))
        return

    outfile.write("# {} {}\n".format(traj_num, len(sputt), ave_surf, rms_val))

    for i in range(len(sputt)):

        x_d = (((sputt[i][0] + box_half_len) % box_len) + x_shift ) % box_len \
              - backward_a


        y_d = (((sputt[i][1] + box_half_len) % box_len) + y_shift ) % box_len \
              - box_half_len

        z_d = sputt[i][2] - ave_surf

        outfile.write("{} {} {}\n".format(x_d, y_d, z_d))




open_str = 'a'
just_test = False
if not just_test:
    rms_file = open("results/rms.txt", open_str)

    sputt_proj_file = open("results/sputtered_proj.txt", open_str)
    sputt_sample_file = open("results/sputtered_sample.txt", open_str)

    # ---
    redist_proj_file = open("results/redist_proj.txt", open_str)
    redist_sum_proj_file = open("results/redist_sum_proj.txt", open_str)

    redist_sample_file = open("results/redist_sample.txt", open_str)
    redist_sum_sample_file = open("results/redist_sum_sample.txt", open_str)


basedir = "/home/dawid/dumps/hobler_redist/si_kr_2kev/dumpse"
dirs = sorted(listdir(basedir))

for dir in dirs[500:]:
    dirr = basedir + "/" + dir
    if just_test:
        print dirr
        break

    traj_num = dir.split(".")[0].split("_")[-1]
    sputt_list, shifted_list, start_z_list, ave_surf, rms_val, projcoord = calc_displacment(dirr)
    #(sputt_list, shifted_list, ave_surf, rms_val, projcoord)

    #print sorted(results[1][0][:,0])[:20]
    #print (sorted(results[1][0][:,0])[-20:])[::-1]
    # 1. write rms results
    rms_file.write("{} {} {}\n".format(traj_num, ave_surf, rms_val))

    # 2. wirte redist sum
    redist_sum_sample_file.write("{} {}\n".format(traj_num, sum(shifted_list[0][:,0])))
    redist_sum_proj_file.write("{} {}\n".format(traj_num, sum(shifted_list[1][:,0])))


    write_redist(traj_num, shifted_list[0], start_z_list[0], redist_sample_file, ave_surf, rms_val)
    write_redist(traj_num, shifted_list[1], start_z_list[1], redist_proj_file, ave_surf, rms_val)

    write_sputt(traj_num, sputt_list[0], sputt_sample_file, ave_surf, rms_val,  projcoord)
    write_sputt(traj_num, sputt_list[1], sputt_proj_file, ave_surf, rms_val,  projcoord)

