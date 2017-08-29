from simpy import analyze
from simpy import write
from scipy.signal import find_peaks_cwt
from scipy.optimize import curve_fit

import os
import numpy as np
import sys

GLOBAL_WRITE = open("/tmp/global.xyz", 'w')

class ResultWriter:
    def __init__(self, input_dir, output_dir, impact_deg, extend=False, basename="", cutoff=0.001):
        self.impact_theta = impact_deg
        self.cutoff=cutoff

        self.redist_interval_dirname = ["{}/{}".format(output_dir, "redist_funz_interval_sample"),
                                        "{}/{}".format(output_dir, "redist_funz_interval_proj")]

        self.displacement_z_interval_dirname = "{}/{}".format(output_dir, "displ_x")
        self.displacement_abs_interval_dirname = "{}/{}".format(output_dir, "displ_abs")

        self.first_moment_histogram_dir = "{}/{}".format(output_dir, "m1_depth_redist")

        self.xhist_sputter_dir = "{}/{}".format(output_dir, "sput_x_histogram")

        self.depth_profile_dirname = "{}/{}".format(output_dir, "atoms_depth")

        if extend:
            write_str = "a"
        else:
            write_str = "w"
        print input_dir
        self.input_dir = input_dir
        self.dirs = sorted(os.listdir(input_dir))

        self.analyze_interval = 150

        # arrays
        self.erosive_ave_arr = [np.zeros((0,3), dtype=float), np.zeros((0,3), dtype=float)]
        self.redist_ave_arr = [np.zeros((0,1), dtype=float), np.zeros((0,1), dtype=float)]

        # dp binning
        lbin_start = -200.0+0.454545454545
        lbin_end = 10.0+0.454545454545
        lbin_size  = 1.35775 #4.0

        lbin_step = int((lbin_end-lbin_start)/lbin_size+1)
        self.dp_bin_centers = np.linspace(lbin_start, lbin_end, lbin_step)
        self.dp_bin_edges = self.dp_bin_centers[:-1]+(0.5*lbin_size)

        # x sputtering binning
        lbin_start = -120.0
        lbin_end = 120.0
        lbin_size  = 1.0

        lbin_step = (lbin_end-lbin_start)/lbin_size+1
        self.xhist_bin_centers = np.linspace(lbin_start, lbin_end, lbin_step)
        self.xhist_bin_edges = self.xhist_bin_centers[:-1]+(0.5*lbin_size)


        # redist funz initialize
        self.bin_start = -100.0
        self.bin_end = 20.0
        self.bin_size = 1.0

        self.bin_steps = (self.bin_end-self.bin_start)/self.bin_size + 1
        self.bin_centers = np.linspace(self.bin_start, self.bin_end, self.bin_steps)
        self.bin_edges = self.bin_centers[:-1]+(0.5*self.bin_size)

        self.bin_total_displacement = [np.zeros(len(self.bin_centers), dtype=float),
                                       np.zeros(len(self.bin_centers), dtype=float)]

        self.bin_local_displacement = [np.zeros(len(self.bin_centers), dtype=float),
                                       np.zeros(len(self.bin_centers), dtype=float)]

        # dislacement cnt
        # z displacement cnt
        self.z_bin_size = 0.1

        self.z_bin_start = -10
        self.z_bin_end = 10

        self.z_bin_steps = (self.z_bin_end-self.z_bin_start)/self.z_bin_size + 1
        self.z_bin_centers = np.linspace(self.z_bin_start, self.z_bin_end, self.z_bin_steps)
        self.z_bin_edges = self.z_bin_centers[:-1]+(0.5*self.z_bin_size)
        self.z_bin_centers = np.power(self.z_bin_centers, 3.0)

        self.z_bin_total_cnt = [np.zeros(len(self.z_bin_centers), dtype=float),
                                np.zeros(len(self.z_bin_centers), dtype=float)]

        self.z_bin_local_cnt = [np.zeros(len(self.z_bin_centers), dtype=float),
                                np.zeros(len(self.z_bin_centers), dtype=float)]

        # first moment histograms
        self.mf_hist_krimplantation = np.zeros(self.bin_centers.shape, dtype=float)
        self.mf_hist_krsput = np.zeros(self.bin_centers.shape, dtype=float)
        self.mf_hist_sisput = np.zeros(self.bin_centers.shape, dtype=float)

        self.mf_hist_krrealoc = np.zeros(self.bin_centers.shape, dtype=float)
        self.mf_hist_sirealoc = np.zeros(self.bin_centers.shape, dtype=float)

        # hist
        self.xhist_kr = np.zeros(self.xhist_bin_centers.shape, dtype=float)
        self.xhist_si = np.zeros(self.xhist_bin_centers.shape, dtype=float)

        # abs displacement cnt
        self.abs_bin_size = 0.1

        self.abs_bin_start = 0
        self.abs_bin_end = 10

        self.abs_bin_steps = (self.abs_bin_end-self.abs_bin_start)/self.abs_bin_size + 1
        self.abs_bin_centers = np.linspace(self.abs_bin_start, self.abs_bin_end, self.abs_bin_steps)
        self.abs_bin_edges = self.abs_bin_centers[:-1]+(0.5*self.abs_bin_size)
        self.abs_bin_centers = np.power(self.abs_bin_centers, 3.0)

        self.abs_bin_total_cnt = [np.zeros(len(self.abs_bin_centers), dtype=float),
                                  np.zeros(len(self.abs_bin_centers), dtype=float)]

        self.abs_bin_local_cnt = [np.zeros(len(self.abs_bin_centers), dtype=float),
                                  np.zeros(len(self.abs_bin_centers), dtype=float)]


        # Kr vs Si volume calculation
        self.ave_surf_dist = 10
        self.volume_height = 20

        # first moment
        self.m_first = open("{}/{}{}".format(output_dir, basename, "first_moment.txt"), write_str)
        if write_str == "w":
            self.m_first.write("traj  kr_impl_mave_{0}  kr_sput_mave_{0} si_sput_mave_{0}  kr_reloc_mave_{0} si_reloc_mave_{0}\n".format(self.analyze_interval))


        self.mf_krimplanted = [np.zeros((0,1), dtype=float), np.zeros((0,1), dtype=float)]
        self.mf_krsputt = [np.zeros((0,1), dtype=float), np.zeros((0,1), dtype=float)]
        self.mf_sisput = [np.zeros((0,1), dtype=float), np.zeros((0,1), dtype=float)]

        self.mf_krrealoc = [np.zeros((0,1), dtype=float), np.zeros((0,1), dtype=float)]
        self.mf_sirealoc = [np.zeros((0,1), dtype=float), np.zeros((0,1), dtype=float)]

        # open files
        self.rms_file = open("{}/{}{}".format(output_dir, basename, "surface.txt"), write_str)
        if write_str == "w":
            self.rms_file.write("traj ave_surf RMS\n")

        self.kr_vs_si = open("{}/{}{}".format(output_dir, basename, "kr_vs_si.txt"), write_str)
        if write_str == "w":
            self.kr_vs_si.write("traj volume volume_z_bottom volume_z_top  kr_cnt si_cnt\n")


        # DP histograms

        '''
        self.m_erosive = [open("{}/{}{}".format(output_dir, basename, "erosive_sample.txt"), write_str),
                          open("{}/{}{}".format(output_dir, basename, "erosive_proj.txt"), write_str)]
        if write_str == "w":
            self.m_erosive[0].write("traj  x y z  x_mave_{0} y_mave_{0} z_mave_{0}\n".format(self.analyze_interval))
            self.m_erosive[1].write("traj  x y z  x_mave_{0} y_mave_{0} z_mave_{0}\n".format(self.analyze_interval))


        self.m_redist_sum = [open("{}/{}{}".format(output_dir, basename, "redistributive_sum_sample.txt"), write_str),
                             open("{}/{}{}".format(output_dir, basename, "redistributive_sum_proj.txt"), write_str)]
        if write_str == "w":
            self.m_redist_sum[0].write("traj M_redist  M_redist_mave_{}\n".format(self.analyze_interval))
            self.m_redist_sum[1].write("traj M_redist  M_residt_mave_{}\n".format(self.analyze_interval))

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


        '''

        self.amorphus_layer = open("{}/{}{}".format(output_dir, basename, "ac_interface.txt"), write_str)
        if write_str == "w":
            self.amorphus_layer.write("traj interface_depth  interface_width  k_param\n")

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
            # projectile is implanted
            ecoord = eframe['coord'][:-1]
            eid = eframe['id'][:-1]
            etype = eframe['type'][:-1]

            self.projectile_implanted = eframe['coord'][-1]
        else:
            # bounced out
            ecoord = eframe['coord']
            eid = eframe['id']
            etype = eframe['type']

            self.projectile_implanted = None

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
        #sframe.remove_type(1)
        #sframe.remove_lesser_than(-35, 2)
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
        self.redist_ave_arr[mode] = np.append(self.redist_ave_arr[mode], sum(self.shifted_list[mode][:,0]))
        out_str = "{}  {}  {}\n".format(self.traj_num, sum(self.shifted_list[mode][:,0]), np.mean(self.redist_ave_arr[mode][self.start_index:]))
        self.m_redist_sum[mode].write(out_str)

    def calc_redist_funz(self, mode):
        data_bin_index = np.digitize(self.start_z_list[mode][:,2]-self.ave_surf, self.bin_edges)
        shifted = self.shifted_list[mode]

        for i in xrange(len(data_bin_index)):
            self.bin_total_displacement[mode][data_bin_index[i]] += shifted[i][0]
            self.bin_local_displacement[mode][data_bin_index[i]] += shifted[i][0]



        if self.traj_num%self.analyze_interval == 0:
            if not os.path.exists(self.redist_interval_dirname[mode]):
                os.makedirs(self.redist_interval_dirname[mode])

            affix_name = "funz_sample.txt" if mode == 0 else "funz_proj.txt"
            start_num = self.traj_num-self.analyze_interval+1
            end_num = self.traj_num
            final_name = "{:04d}-{:04d}_{}".format(start_num, end_num, affix_name)
            local_redist_file = open("{}/{}".format(self.redist_interval_dirname[mode], final_name), 'w')
            local_redist_file.write("dz_from_ave_surf M_redist\n")


            for center, displ in zip(self.bin_centers, self.bin_local_displacement[mode]/self.analyze_interval):
                local_redist_file.write("{} {}\n".format(center, displ))

            self.bin_local_displacement[mode].fill(0.0)


    def write_redist_funz(self, mode):
        for center, displ in zip(self.bin_centers, self.bin_total_displacement[mode]):
            self.m_redist_funz[mode].write("{} {}\n".format(center, displ/self.analyzed_total))


    def calc_displ(self):
        # z
        # si
        squbed = (abs(self.shifted_list[0][:,0])/self.shifted_list[0][:,0])*np.power(np.abs(self.shifted_list[0][:,0]), 1./3.)
        data_z_bin_index = np.digitize(squbed, self.z_bin_edges)

        for bin_index in data_z_bin_index:
            self.z_bin_local_cnt[0][bin_index] += 1

        # kr
        squbed = (abs(self.shifted_list[1][:,0])/self.shifted_list[1][:,0])*np.power(np.abs(self.shifted_list[1][:,0]), 1./3.)
        data_z_bin_index = np.digitize(squbed, self.z_bin_edges)

        for bin_index in data_z_bin_index:
            self.z_bin_local_cnt[1][bin_index] += 1


        # abs
        # si
        abs_displ = np.power(np.linalg.norm(self.shifted_list[0],axis=1), 1.0/3.0)
        data_abs_bin_index = np.digitize(abs_displ, self.abs_bin_edges)

        for bin_index in data_abs_bin_index:
            self.abs_bin_local_cnt[0][bin_index] += 1

        # kr
        abs_displ = np.power(np.linalg.norm(self.shifted_list[1],axis=1), 1.0/3.0)
        data_abs_bin_index = np.digitize(abs_displ, self.abs_bin_edges)

        for bin_index in data_abs_bin_index:
            self.abs_bin_local_cnt[1][bin_index] += 1


        if self.traj_num%self.analyze_interval == 0:
            if not os.path.exists(self.displacement_z_interval_dirname):
                os.makedirs(self.displacement_z_interval_dirname)


            if not os.path.exists(self.displacement_abs_interval_dirname):
                os.makedirs(self.displacement_abs_interval_dirname)

            start_num = self.traj_num-self.analyze_interval+1
            end_num = self.traj_num

            final_z_name = "{:04d}-{:04d}.txt".format(start_num, end_num)
            final_abs_name = "{:04d}-{:04d}.txt".format(start_num, end_num)

            local_z_file = open("{}/{}".format(self.displacement_z_interval_dirname, final_z_name),"w")
            local_z_file.write("x_displ n_Kr n_Si\n")
            local_abs_file = open("{}/{}".format(self.displacement_abs_interval_dirname, final_abs_name),"w")
            local_abs_file.write("abs_displ n_Kr n_Si\n")

            for center, cnt_kr, cnt_si in zip(self.z_bin_centers, self.z_bin_local_cnt[1], self.z_bin_local_cnt[0]):
                local_z_file.write("{} {} {}\n".format(center, cnt_kr/self.analyze_interval, cnt_si/self.analyze_interval))

            for center, cnt_kr, cnt_si in zip(self.abs_bin_centers, self.abs_bin_local_cnt[1], self.abs_bin_local_cnt[0]):
                local_abs_file.write("{} {} {}\n".format(center, cnt_kr/self.analyze_interval, cnt_si/self.analyze_interval))

            self.z_bin_local_cnt[0].fill(0.0)
            self.z_bin_local_cnt[1].fill(0.0)
            self.abs_bin_local_cnt[0].fill(0.0)
            self.abs_bin_local_cnt[1].fill(0.0)

    def write_displ(self, mode):
        # z
        for center, cnt in zip(self.z_bin_centers, self.z_bin_total_cnt[mode]):
            self.z_displ[mode].write("{} {}\n".format(center, cnt/self.analyzed_total))

        # abs
        for center, cnt in zip(self.abs_bin_centers, self.abs_bin_total_cnt[mode]):
            self.abs_displ[mode].write("{} {}\n".format(center, cnt/self.analyzed_total))

    def write_first_moment(self):
        # compute projectile hit place and shift requirements
        box_len_x = -self.box[0][0] + self.box[0][1]
        box_len_y = -self.box[1][0] + self.box[1][1]

        box_half_len_x = box_len_x*0.5
        box_half_len_y = box_len_y*0.5

        backward_x = 0.25 * box_len_x

        # computing projectile x path - from start to average sureface crossing
        z_proj_path = self.projcoord[2] - self.ave_surf
        x_proj_path = z_proj_path * np.tan(np.radians(self.impact_theta))

        # moving xlow boundary of sample to 0 and calcualting impact point wrapping
        x_hit_position = ((x_proj_path + self.projcoord[0]) - self.box[0][0]) % box_len_x
        x_hit_position_unwrapped = x_proj_path + self.projcoord[0]
        # calculating how much we have to move sample to get hit position at 0.25*len_of_the_sample
        x_shift = backward_x - x_hit_position

        # how much we have to move y in oreder to be in center
        y_hit_position = (self.projcoord[1] - self.box[1][0]) % box_len_y
        y_shift = box_half_len_y - y_hit_position

        # Kr implanted (x) // x1-x0 self.projectile_implanted

        # !!!!
        if self.projectile_implanted is None:
            self.mf_krimplanted = np.append(self.mf_krimplanted, 0.0)
        else:
            implanted_distance = self.projectile_implanted[0]-x_hit_position_unwrapped
            self.mf_krimplanted = np.append(self.mf_krimplanted, implanted_distance)
            proj_indexes = np.digitize(self.projectile_implanted[2] - self.ave_surf, self.bin_edges)
            self.mf_hist_krimplantation[proj_indexes] += implanted_distance

        # 0 - Si  1 - Kr
        # Kr sputtered (x) // x0-x1
        if self.sputt_list[1] is None:
            self.mf_krsputt = np.append(self.mf_krsputt, 0.0)
        else:
            kr_sputtered_total = 0.0
            for x_start_pos in self.sputt_list[1][:,0]:
                kr_sputtered_total -= (((x_start_pos- self.box[0][0]) % box_len_x) + x_shift ) % box_len_x - backward_x

            self.mf_krsputt = np.append(self.mf_krsputt, kr_sputtered_total)

            start_pos = self.sputt_list[1]
            proj_indexes = np.digitize(start_pos[:,2] - self.ave_surf, self.bin_edges)
            x_displ_arr = []
            for i in xrange(len(proj_indexes)):
                displ_value = (((start_pos[i][0]- self.box[0][0]) % box_len_x) + x_shift ) % box_len_x - backward_x
                self.mf_hist_krsput[proj_indexes[i]] -= displ_value #(((start_pos[i][0]- self.box[0][0]) % box_len_x) + x_shift ) % box_len_x - backward_x

                x_displ_arr.append(displ_value)


            proj_indexes = np.digitize(np.array(x_displ_arr), self.xhist_bin_edges)
            for indexes in proj_indexes:
                self.xhist_kr[indexes] += 1



        # Si sputtered (x) // x0-x1
        if self.sputt_list[0] is None:
            self.mf_sisput = np.append(self.mf_sisput, 0.0)
        else:
            si_sputtered_total = 0.0
            #GLOBAL_WRITE.write("{}\n{}\n".format(len(self.sputt_list[0]), "name"))
            for atomi in range(len(self.sputt_list[0])):
                x_sputtered_pos = (((self.sputt_list[0][atomi][0] - self.box[0][0]) % box_len_x) + x_shift ) % box_len_x - backward_x

                y_sputtered_pos = (((self.sputt_list[0][atomi][1] - self.box[1][0]) % box_len_y) + y_shift ) % box_len_y - box_half_len_y

                z_sputtered_pos = self.sputt_list[0][atomi][2] - self.ave_surf

                GLOBAL_WRITE.write("Si {} {} {}\n".format(x_sputtered_pos, y_sputtered_pos, z_sputtered_pos))
                si_sputtered_total -= x_sputtered_pos

            self.mf_sisput = np.append(self.mf_sisput, si_sputtered_total)

            start_pos = self.sputt_list[0]
            proj_indexes = np.digitize(start_pos[:,2] - self.ave_surf, self.bin_edges)
            x_displ_arr = []
            for i in xrange(len(proj_indexes)):
                displ_value = (((start_pos[i][0]- self.box[0][0]) % box_len_x) + x_shift ) % box_len_x - backward_x
                self.mf_hist_sisput[proj_indexes[i]] -= displ_value

                x_displ_arr.append(displ_value)

            proj_indexes = np.digitize(np.array(x_displ_arr), self.xhist_bin_edges)
            for indexes in proj_indexes:
                self.xhist_si[indexes] += 1



        # Kr relocated  // x2-x1
        if self.shifted_list[1]  is None:
            lsum = 0.0
        else:
            lsum = sum(self.shifted_list[1][:,0])

            proj_indexes = np.digitize(self.start_z_list[1][:,2] - self.ave_surf, self.bin_edges)
            shift_list = self.shifted_list[1][:,0]
            for i in xrange(len(proj_indexes)):
                self.mf_hist_krrealoc[proj_indexes[i]] += shift_list[i]

        self.mf_krrealoc = np.append(self.mf_krrealoc, lsum)

        # Si relocated // x2-x1
        if self.shifted_list[0]  is None:
            lsum = 0.0
        else:
            lsum = sum(self.shifted_list[0][:,0])

            proj_indexes = np.digitize(self.start_z_list[0][:,2] - self.ave_surf, self.bin_edges)
            shift_list = self.shifted_list[0][:,0]
            for i in xrange(len(proj_indexes)):
                self.mf_hist_sirealoc[proj_indexes[i]] += shift_list[i]
        self.mf_sirealoc = np.append(self.mf_sirealoc, lsum)

        # write results
        rs = max(self.analyzed_total-1-self.analyze_interval, 0)
        re = self.analyzed_total+1

        krimpl = np.average(self.mf_krimplanted[rs:re])
        krsput = np.average(self.mf_krsputt[rs:re])
        sisput = np.average(self.mf_sisput[rs:re])
        krre = np.average(self.mf_krrealoc[rs:re])
        sire = np.average(self.mf_sirealoc[rs:re])

        w_stirng = "{} {}  {} {}  {} {}\n".format(self.traj_num, krimpl,  krsput, sisput, krre, sire)

        self.m_first.write(w_stirng)

        # /////
        if self.traj_num%self.analyze_interval == 0:
            if not os.path.exists(self.first_moment_histogram_dir):
                os.makedirs(self.first_moment_histogram_dir)

            start_num = self.traj_num-self.analyze_interval+1
            end_num = self.traj_num
            local_hist_name = "{:04}-{:04}.txt".format(start_num, end_num)

            local_hist_file =  open("{}/{}".format(self.first_moment_histogram_dir, local_hist_name), 'w')
            local_hist_file.write("depth-av_surf kr_impl kr_sput si_sput  kr_reloc si_reloc\n")

            div = self.analyze_interval
            for i in range(len(self.bin_centers)):
                outstr = "{} {}  {} {}  {} {}\n".format(self.bin_centers[i], self.mf_hist_krimplantation[i]/div,
                                                 self.mf_hist_krsput[i]/div, self.mf_hist_sisput[i]/div,
                                                 self.mf_hist_krrealoc[i]/div, self.mf_hist_sirealoc[i]/div)
                local_hist_file.write(outstr)



            self.mf_hist_krimplantation.fill(0.0)

            self.mf_hist_krsput.fill(0.0)
            self.mf_hist_sisput.fill(0.0)

            self.mf_hist_krrealoc.fill(0.0)
            self.mf_hist_sirealoc.fill(0.0)


            if not os.path.exists(self.xhist_sputter_dir):
                os.makedirs(self.xhist_sputter_dir)

            local_xhist_file =  open("{}/{}".format(self.xhist_sputter_dir, local_hist_name), 'w')
            local_xhist_file.write("x2-x1  n_kr  n_si\n")


            for i in range(len(self.xhist_bin_centers)):
                outstr = "{} {} {}\n".format(self.xhist_bin_centers[i], self.xhist_kr[i]/div, self.xhist_si[i]/div)
                local_xhist_file.write(outstr)

            self.xhist_si.fill(0.0)
            self.xhist_kr.fill(0.0)


    def write_depth_profile(self):

        if self.traj_num%self.analyze_interval != 1:
            return

        if not os.path.exists(self.depth_profile_dirname):
            os.makedirs(self.depth_profile_dirname)

        final_name = "{}.txt".format(self.traj_num-1)
        local_file = open("{}/{}".format(self.depth_profile_dirname, final_name), 'w')
        local_file.write("depth-ave_surf n_Kr n_Si\n")

        indexed_si = np.digitize(self.start_all_z_list[0]- self.ave_surf, self.dp_bin_edges)
        binned_si = np.zeros(self.dp_bin_centers.shape, dtype=int)
        for i in indexed_si:
            binned_si[i] += 1

        indexed_kr = np.digitize(self.start_all_z_list[1]- self.ave_surf, self.dp_bin_edges)
        binned_kr = np.zeros(self.dp_bin_centers.shape, dtype=int)
        for i in indexed_kr:
            binned_kr[i] += 1

        for center, kr, si in zip(self.dp_bin_centers, binned_kr, binned_si):
            local_file.write("{}  {} {}\n".format(center, kr, si))

    def write_erosive(self, mode):
        box_len_x = -self.box[0][0] + self.box[0][1]
        box_len_y = -self.box[1][0] + self.box[1][1]

        box_half_len_x = box_len_x*0.5
        box_half_len_y = box_len_y*0.5

        backward_x = 0.25 * box_len_x

        z_proj_path = self.projcoord[2] - self.ave_surf
        x_proj_path = z_proj_path * np.tan(np.radias(self.impact_theta))

        x_hit_position = ((x_proj_path + self.projcoord[0]) - self.box[0][0]) % box_len_x
        x_shift = backward_x - x_hit_position

        y_hit_position = (self.projcoord[1] - self.box[1][0]) % box_len_y
        y_shift = box_half_len_y - y_hit_position


        if self.sputt_list[mode] is None:
            self.erosive_ave_arr[mode] =  np.append(self.erosive_ave_arr[mode], np.array([[0.,0.,0.]]), axis=0)
            print "skipping sputtering - None"
            out_str = "{}  {} {} {}".format(self.traj_num, 0.0, 0.0, 0.0)
            out_str = "{}  {} {} {}\n".format(out_str, *np.mean(self.erosive_ave_arr[mode][self.start_index:], axis=0))
            self.m_erosive[mode].write(out_str)
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
        self.erosive_ave_arr[mode] = np.append(self.erosive_ave_arr[mode], np.array([[x_tot, y_tot, z_tot]]), axis=0)
        out_str = "{}  {} {} {}".format(self.traj_num, x_tot, y_tot, z_tot)
        out_str = "{}  {} {} {}\n".format(out_str, *np.mean(self.erosive_ave_arr[mode][self.start_index:], axis=0))
        self.m_erosive[mode].write(out_str)

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
        try:
            coef, pcov  = curve_fit(sigmoid, xdata, ydata, p0=init_optimal)
        except:
            coef = [0, np.inf, 0,0]

        fit_curve = sigmoid(xdata, coef[0], coef[1], coef[2], coef[3])

        if False:
            tmp_write = np.transpose(np.array([xdata, ydata, fit_curve]))
            np.savetxt("/tmp/curve_fit.txt", tmp_write)
            np.savetxt("/tmp/hist_data.txt", np.transpose(np.array([bin_centers, bin_cnt])))

        amorph_width = self.ave_surf - coef[0]
        amorph_interface_width = (-1.0/coef[1])*np.log((1-cutoff_value)/cutoff_value)

        self.amorphus_layer.write("{}  {} {} {}\n".format(self.traj_num, -amorph_width, amorph_interface_width, coef[1]))

    def run(self, runrange=False, skip=1):
        if not runrange:
            runrange = [0,100000]

        self.analyzed_total = 0
        for crdir in self.dirs[runrange[0]:runrange[1]:skip]:
            dirr = "{}/{}".format(self.input_dir, crdir)
            self.traj_num = int(crdir.split(".")[0].split("_")[-1])

            self.strat_index = 1 + self.analyzed_total - self.analyze_interval
            if self.strat_index < 0:
                self.start_index = 0


            self.calc_displacment(dirr)

            self.rms_file.write("{} {} {}\n".format(self.traj_num, self.ave_surf, self.rms_val))
            self.write_amorphous_layer()
            self.write_first_moment()
            self.write_depth_profile()

            self.calc_displ()

            self.write_kr_vs_si()

            #self.write_displ()

            '''
            self.write_erosive(0)
            self.write_erosive(1)

            self.write_redist_sum(0)
            self.write_redist_sum(1)

            # TO DO: write periodic save
            # DONE ???
            self.calc_redist_funz(0)
            self.calc_redist_funz(1)

            # TO DO: write periodic save
            self.calc_displ(0)
            self.calc_displ(1)


            self.write_amorphous_layer()
            '''

            self.analyzed_total += 1

        '''
        self.write_redist_funz(0)
        self.write_redist_funz(1)

        self.write_displ(0)
        self.write_displ(1)
        '''

reswrite = ResultWriter("/home/dawid/dumps/hobler_redist/diff_deg/dumpse_40",
                        "/tmp/test",
                        40)
#reswrite.run([899, 1201])
reswrite.run([650, 1000])
print "test"
#reswrite.run([1501, 5000])

