from simpy import analyze
from simpy import write
from scipy.signal import find_peaks_cwt
from scipy.optimize import curve_fit

from os import listdir
import numpy as np
import sys

class ResultWriter:
    def __init__(self, input_dir, output_dir, extend=False, basename="", cutoff=0.001):
        self.cutoff=cutoff

        if extend:
            write_str = "a"
        else:
            write_str = "w"
        print input_dir
        self.input_dir = input_dir
        self.dirs = sorted(listdir(input_dir))

        self.analyzed_total = 0

        # redist funz initialize
        self.bin_start = -80.0
        self.bin_end = 10.0
        self.bin_size = 1.0

        self.bin_steps = (self.bin_end-self.bin_start)/self.bin_size + 1
        self.bin_centers = np.linspace(self.bin_start, self.bin_end, self.bin_steps)
        self.bin_edges = self.bin_centers[:-1]+(0.5*self.bin_size)

        self.bin_total_displacement = [np.zeros(len(self.bin_centers), dtype=float),
                                       np.zeros(len(self.bin_centers), dtype=float)]

        # dislacement cnt
        # z displacement cnt
        self.z_bin_size = 1

        self.z_bin_start = -80
        self.z_bin_end = 80

        self.z_bin_steps = (self.z_bin_end-self.z_bin_start)/self.z_bin_size + 1
        self.z_bin_centers = np.linspace(self.z_bin_start, self.z_bin_end, self.z_bin_steps)
        self.z_bin_edges = self.z_bin_centers[:-1]+(0.5*self.z_bin_size)

        self.z_bin_total_cnt = [np.zeros(len(self.z_bin_centers), dtype=float),
                                np.zeros(len(self.z_bin_centers), dtype=float)]

        # abs displacement cnt
        self.abs_bin_size = 1

        self.abs_bin_start = 0
        self.abs_bin_end = 200

        self.abs_bin_steps = (self.abs_bin_end-self.abs_bin_start)/self.abs_bin_size + 1
        self.abs_bin_centers = np.linspace(self.abs_bin_start, self.abs_bin_end, self.abs_bin_steps)
        self.abs_bin_edges = self.abs_bin_centers[:-1]+(0.5*self.abs_bin_size)

        self.abs_bin_total_cnt = [np.zeros(len(self.abs_bin_centers), dtype=float),
                                  np.zeros(len(self.abs_bin_centers), dtype=float)]


        # Kr vs Si volume calculation
        self.ave_surf_dist = 10
        self.volume_height = 20

        # open files
        self.rms_file = open("{}/{}{}".format(output_dir, basename, "rms.txt"), write_str)
        if write_str == "w":
            self.rms_file.write("traj ave_surf RMS\n")

        self.m_erosive = [open("{}/{}{}".format(output_dir, basename, "erosive_sample.txt"), write_str),
                          open("{}/{}{}".format(output_dir, basename, "erosive_proj.txt"), write_str)]
        if write_str == "w":
            self.m_erosive[0].write("traj x y z\n")
            self.m_erosive[1].write("traj x y z\n")


        self.m_redist_sum = [open("{}/{}{}".format(output_dir, basename, "redistributive_sum_sample.txt"), write_str),
                             open("{}/{}{}".format(output_dir, basename, "redistributive_sum_proj.txt"), write_str)]
        if write_str == "w":
            self.m_redist_sum[0].write("traj M_redist\n")
            self.m_redist_sum[1].write("traj M_redist\n")

        self.m_redist_funz = [open("{}/{}{}".format(output_dir, basename, "redistributive_funz_sample.txt"), write_str),
                              open("{}/{}{}".format(output_dir, basename, "redistributive_funz_proj.txt"), write_str)]
        if write_str == "w":
            self.m_redist_funz[0].write("dz_from_ave_surf M_redist\n")
            self.m_redist_funz[1].write("dz_from_ave_surf M_redist\n")

        self.z_displ = [open("{}/{}{}".format(output_dir, basename, "displacement_z_sample.txt"), write_str),
                        open("{}/{}{}".format(output_dir, basename, "displacement_z_proj.txt"), write_str)]
        if write_str == "w":
            self.z_displ[0].write("z_displacement atom_cnt\n")
            self.z_displ[1].write("z_displacement atom_cnt\n")

        self.abs_displ =[open("{}/{}{}".format(output_dir, basename, "displacement_abs_sample.txt"), write_str),
                         open("{}/{}{}".format(output_dir, basename, "displacement_abs_proj.txt"), write_str)]
        if write_str == "w":
            self.abs_displ[0].write("abs_displacement atom_cnt\n")
            self.abs_displ[1].write("abs_displacement atom_cnt\n")

        self.kr_vs_si = open("{}/{}{}".format(output_dir, basename, "kr_vs_si.txt"), write_str)
        if write_str == "w":
            self.kr_vs_si.write("traj volume volume_z_bottom volume_z_top  kr_cnt si_cnt\n")

        self.amorphus_layer = open("{}/{}{}".format(output_dir, basename, "amorphous_layer_param.txt"), write_str)
        if write_str == "w":
            self.amorphus_layer.write("traj layer_width  interface_width  k_param\n")

    def calc_displacment(self, dir):
        traj = analyze.Traj(dir)

        # remove new ion form start and from end if implanted
        sframe = traj.read()
        scoord = sframe['coord'][:-1]
        projectile_id = sframe['id'][-1]
        sid = sframe['id'][:-1]
        stype = sframe['type'][:-1]

        self.amorph_coord = scoord[:, 2]
        self.start_all_z_list = [scoord[stype==2][:,2],
                                 scoord[stype==1][:,2]]
        self.box = sframe['box']


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

        #start_z = np.empty(len(ecoord), dtype=float)
        start_z = np.empty((len(ecoord), 3), dtype=float)

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

            #start_z[iend] = scoord[i][2]
            start_z[iend] = scoord[i]

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
        cutoff = self.cutoff
        above_cutoff = np.abs(shift[:,0]) >= cutoff

        shifted_list = (shift[np.logical_and(above_cutoff, shift_si_flag)],
                  shift[np.logical_and(above_cutoff, np.logical_not(shift_si_flag))])

        start_z_list = (start_z[np.logical_and(above_cutoff, shift_si_flag)],
                  start_z[np.logical_and(above_cutoff, np.logical_not(shift_si_flag))])

        if iend != len(ecoord):
            print "ERROR !!!!"
            print "missing atoms"
            print dir
            sys.exit()

        # ------------------------
        # RMS and surf lvl analyze
        # ------------------------
        sframe.remove_type(1)
        sframe.remove_lesser_than(-35, 2)
        old_box = sframe['box']

        sframe.periodic([[-1,1], [-1,1], [0,0]])
        rms = analyze.RMS(sframe['coord'], 2.1, 2.1, old_box)

        ave_surf, rms_val = rms.compute()

        self.sputt_list = sputt_list
        self.shifted_list = shifted_list
        self.start_z_list = start_z_list
        self.ave_surf = ave_surf
        self.rms_val = rms_val
        self.projcoord = projcoord
        return
        #return (sputt_list, shifted_list, start_z_list, ave_surf, rms_val, projcoord)

        #bin_size = 2
        #for i in range(-18,5)[::-1]:
        #    start = ave_surf + i*bin_size
        #    stop = start - bin_size
        #    layer = np.logical_and(start_z <= start, start_z > stop)
        #    print "l ({}, {}) : {}   // atoms = {}".format(start, stop, sum(shift[layer][:,2]), sum(layer))
        #name = dir.split('/')[-1].split(".")[0] + ".xyz"
        #write.xyz(name).write(rms.get_frame())


    def get_rms(self, frame):
        frame.remove_type(1)
        frame.remove_lesser_than(-30,2)
        old_box = frame['box']

        frame.periodic([[-1,1], [-1,1], [0,0]])
        rms = analyze.RMS(frame['coord'], 2.1, 2.1, old_box)

        rms.compute()
        #frame = rms.get_frame()
        #write.xyz("out.xyz").write(frame)


    #def write_redist(traj_num, redist, start_z_list, outfile, ave_surf, rms_val, projcoord=None):
    def write_redist_sum(self, mode):
        self.m_redist_sum[mode].write("{} {}\n".format(self.traj_num, sum(self.shifted_list[mode][:,0])))

    def calc_redist_funz(self, mode):
        data_bin_index = np.digitize(self.start_z_list[mode][:,2]-self.ave_surf, self.bin_edges)
        shifted = self.shifted_list[mode]

        for i in xrange(len(data_bin_index)):
            self.bin_total_displacement[mode][data_bin_index[i]] += shifted[i][0]

    def write_redist_funz(self, mode):
        for center, displ in zip(self.bin_centers, self.bin_total_displacement[mode]):
            self.m_redist_funz[mode].write("{} {}\n".format(center, displ/self.analyzed_total))


    def calc_displ(self, mode):
        # z
        data_z_bin_index = np.digitize(self.shifted_list[mode][:,0], self.z_bin_edges)

        for bin_index in data_z_bin_index:
            self.z_bin_total_cnt[mode][bin_index] += 1

        # abs
        data_abs_bin_index = np.digitize(np.linalg.norm(self.shifted_list[mode],axis=1), self.abs_bin_edges)

        for bin_index in data_abs_bin_index:
            self.abs_bin_total_cnt[mode][bin_index] += 1

    def write_displ(self, mode):
        # z
        for center, cnt in zip(self.z_bin_centers, self.z_bin_total_cnt[mode]):
            self.z_displ[mode].write("{} {}\n".format(center, cnt/self.analyzed_total))

        # abs
        for center, cnt in zip(self.abs_bin_centers, self.abs_bin_total_cnt[mode]):
            self.abs_displ[mode].write("{} {}\n".format(center, cnt/self.analyzed_total))


    def write_erosive(self, mode):
        box_len_x = -self.box[0][0] + self.box[0][1]
        box_len_y = -self.box[1][0] + self.box[1][1]

        box_half_len_x = box_len_x*0.5
        box_half_len_y = box_len_y*0.5

        backward_x = 0.25 * box_len_x

        z_proj_path = self.projcoord[2] - self.ave_surf
        x_proj_path = np.sqrt(3) * z_proj_path

        x_hit_position = ((x_proj_path + self.projcoord[0]) - self.box[0][0]) % box_len_x
        x_shift = backward_x - x_hit_position

        y_hit_position = (self.projcoord[1] - self.box[1][0]) % box_len_y
        y_shift = box_half_len_y - y_hit_position

        if self.sputt_list[mode] is None:
            print "skipping sputtering - None"
            self.m_erosive[mode].write("{}  {} {} {}\n".format(self.traj_num, 0.0, 0.0, 0.0))
            return

        #outfile = self.m_erosive[mode]
        #outfile.write("# {} {}\n".format(self.traj_num, len(self.sputt_list[mode]), self.ave_surf, self.rms_val))

        sputt = self.sputt_list[mode]

        x_tot = y_tot = z_tot = 0.0
        for i in range(len(sputt)):

            x_d = (((sputt[i][0] - self.box[0][0]) % box_len_x) + x_shift ) % box_len_x \
                  - backward_x
            x_tot += x_d


            y_d = (((sputt[i][1] - self.box[1][0]) % box_len_y) + y_shift ) % box_len_y \
                  - box_half_len_y
            y_tot += y_d

            z_d = sputt[i][2] - self.ave_surf
            z_tot += z_d

            #outfile.write("{} {} {}\n".format(x_d, y_d, z_d))
        self.m_erosive[mode].write("{}  {} {} {}\n".format(self.traj_num, x_tot, y_tot, z_tot))

    def write_kr_vs_si(self):
        volume = (self.box[0][1]-self.box[0][0])*(self.box[1][1]-self.box[1][0])*self.volume_height
        volume_top = self.ave_surf - self.ave_surf_dist
        volume_bottom = volume_top - self.volume_height

        si_cnt = sum(np.logical_and(self.start_all_z_list[0]<volume_top, self.start_all_z_list[0]>volume_bottom))
        kr_cnt = sum(np.logical_and(self.start_all_z_list[1]<volume_top, self.start_all_z_list[1]>volume_bottom))

        self.kr_vs_si.write("{} {} {} {}  {} {}\n".format(self.traj_num, volume, volume_bottom, volume_top, kr_cnt, si_cnt))

    def write_amorphous_layer(self):

        def sigmoid(x, x0, k, a, b):
            return  a / (1.0 + np.exp(-k * (x - x0))) + b

        cutoff_value = 0.2

        lo = -160
        hi = 10
        bin_size = 0.1

        bin_number = (hi - lo) / float(bin_size) + 1

        bin_centers = np.linspace(lo, hi, bin_number)
        bin_edges = bin_centers[:-1] + 0.5 * bin_size

        bin_indexes = np.digitize(self.amorph_coord, bin_edges)

        bin_cnt = np.zeros(bin_centers.shape, dtype=int)
        for i in bin_indexes:
            bin_cnt[i] += 1

        # smooth peaks
        bin_cnt = (bin_cnt + np.roll(bin_cnt, 1) + np.roll(bin_cnt, -1))/3.0

        # remove last one for safety
        peaks = find_peaks_cwt(bin_cnt, np.arange(5,10))[:-1]

        xdata = bin_centers[peaks]
        ydata = bin_cnt[peaks]

        init_optimal = [-30, -0.1, 300, 100]
        coef, pcov  = curve_fit(sigmoid, xdata, ydata, p0=init_optimal)

        fit_curve = sigmoid(xdata, coef[0], coef[1], coef[2], coef[3])

        if True:
            tmp_write = np.transpose(np.array([xdata, ydata, fit_curve]))
            np.savetxt("/tmp/curve_fit.txt", tmp_write)
            np.savetxt("/tmp/hist_data.txt", np.transpose(np.array([bin_centers, bin_cnt])))


        amorph_width = self.ave_surf - coef[0]
        amorph_interface_width = (-1.0/coef[1])*np.log((1-cutoff_value)/cutoff_value)

        self.amorphus_layer.write("{}  {} {} {}\n".format(self.traj_num, amorph_width, amorph_interface_width, coef[1]))

    def run(self, runrange=False, skip=1):
        if not runrange:
            runrange = [0,100000]

        for crdir in self.dirs[runrange[0]:runrange[1]:skip]:
            dirr = "{}/{}".format(self.input_dir, crdir)
            self.traj_num = crdir.split(".")[0].split("_")[-1]

            self.calc_displacment(dirr)

            self.rms_file.write("{} {} {}\n".format(self.traj_num, self.ave_surf, self.rms_val))

            # TO DO: write average
            self.write_erosive(0)
            self.write_erosive(1)

            # TO DO: write average
            self.write_redist_sum(0)
            self.write_redist_sum(1)

            # TO DO: write average better redistrubution
            self.calc_redist_funz(0)
            self.calc_redist_funz(1)

            # TO DO: write average
            # better bining methond (x^3)
            self.calc_displ(0)
            self.calc_displ(1)

            self.write_kr_vs_si()

            self.write_amorphous_layer()

            self.analyzed_total += 1

        self.write_redist_funz(0)
        self.write_redist_funz(1)

        self.write_displ(0)
        self.write_displ(1)

reswrite = ResultWriter("/home/dawid/dumps/hobler_redist/diff_deg/dumpse_80", "/tmp")
reswrite.run([750,751], 1)

