"""
This module is used to analyze simulation results.
"""
from simpy import tools
import numpy
import math
import logging
log = logging.getLogger(__name__)

def get(frame, field):
    """
    Fetch array corresponding to require field.
    For singe frame returns fingle array.
    For list of frames returns list of array.

    Args:
        frame (list): list of frames or single frame
        field (str): name of field to fetch

    Returns:
        numpy.ndarray: if numeric field (coords, velocity ...)
        list [str]: if text field (elements)
    """

    try:
        if type(frame) == list and type(frame[0]) == dict:
            fields = frame[0]['format']
            index = frame[0]['format'].index(field) + 1
            return frame[index]
        elif type(frame) == list and type(frame[0]) == list and type(frame[0][0]) == dict:
            fields = frame[0][0]['format']
            out_array = []
            for i in xrange(len(frame)):
                index = frame[i][0]['format'].index(field) + 1
                out_array.append(frame[i][index])
            return out_array
        else:
            log.warning("This is not a frame or a list of frames."
                        "This function accepts returns from analyze.Traj.return()")
    except Exception as exp:
        print type(exp)
        print exp
        log.warning("Something went wong when getting field: %s" % field)
        log.warning("Maybe field does not exist. Available fields:\n%s"
                    % str(fields))

    return None


class Traj:
    """
    Tool for reading simulation dumps with different extensions

    Args:
        filename (str): file name
        format (str): file format
    """

    def __init__(self, filename, format='lammpstrj'):
        try:
            self.infile = open(filename, 'r')
            log.info("Successfully opened file: %s" % filename)
        except IOError:
            log.error("Could not find/read file: %s" % filename)

        function_name = "_read_" + format
        try:
            self.call_read = getattr(self, function_name)
        except AttributeError:
            log.error("Unknow file format: %s" % str(format))

    def read(self, read_n=None):
        """
        Reads trajectory file frame by frame

        Args:
            read_n (int): if positive: number of trajestories to read
                if negative: read the whole file
        Returns:
            frame  if read_n=None
            list of frames  if reand_n=int
        """

        if read_n is None:
            log.info("Rading single frame")
            return self.call_read(False)

        read_whole = False if read_n > 0 else True
        frame = True

        frames = []
        log.info("Read %i frames" % read_n)
        frames_readed = 0
        read_begining = str(read_n) if read_n > 0 else "END"
        while frame and (read_n or read_whole):
            frame = self.call_read(False)
            if frame:
                frames_readed += 1
                frames.append(frame)
                log.info("Current frame (%i / %s)" %
                         (frames_readed, read_begining))

            read_n -= 1

        return frames

    def skip(self, skip_n=1):
        not_end = 1
        skip_begining = skip_n
        log.info("Skip %i frames" % skip_n)
        while not_end and skip_n:
            not_end = self.call_read(True)
            skip_n -= 1

        if not not_end:
            skipped = skip_begining - skip_n
            log.warning("After %i frames skipped end of file encountered" %
                        skipped)

    def _read_lammpstrj(self, skip):
        """
        Args:
            skip (bool): True - skip one frame, False - read and return frame

        Returns:
            frame: single lammps frame
        """

        infile = self.infile
        frame_info = {}

        line = infile.readline()
        if not line:
            return False

        # loop over all header lines
        key_field_name = line.split()[1]
        while key_field_name != "ATOMS":
            if key_field_name == "TIMESTEP":
                timestep = int(infile.readline().split()[0])
                frame_info['timestep'] = int(timestep)

            elif key_field_name == "NUMBER":
                frame_atoms = int(infile.readline().split()[0])
                frame_info['atoms'] = int(frame_atoms)

            elif key_field_name == "BOX":
                box_size = []
                sl = infile.readline().split()
                box_size.append([float(sl[0]), float(sl[1])])
                sl = infile.readline().split()
                box_size.append([float(sl[0]), float(sl[1])])
                sl = infile.readline().split()
                box_size.append([float(sl[0]), float(sl[1])])

                frame_info['box'] = box_size

            elif key_field_name == "TIME":
                time_real = float(infile.readline().split()[0])
                frame_info['time'] = time_real

            else:
                log.warn("Unknow key field name: %s" % key_field_name)

            line = infile.readline()
            key_field_name = line.split()[1]

        # Prepare memory strusture on the header basis
        sp = line.split()
        key_fields = sp[2:]

        # key field groups (coords and velocity)
        # flags value meaning
        #   False - not set
        #   1 - set but memory unallocated
        #   2 - set and memory is allocated
        if 'x' in key_fields and 'y' in key_fields and 'z' in key_fields:
            coords_group = 1
        else:
            coords_group = False

        if 'xs' in key_fields and 'ys' in key_fields and 'zs' in key_fields:
            coords_s_group = 1
        else:
            coords_s_group = False

        if 'xu' in key_fields and 'yu' in key_fields and 'zu' in key_fields:
            coords_u_group = 1
        else:
            coords_u_group = False

        if 'vx' in key_fields and 'vy' in key_fields and 'vz' in key_fields:
            velocity_group = 1
        else:
            velocity_group = False

        frame_info['format'] = []

        # transform key position index in file to index in output list
        key_proxy = [None]*len(key_fields)

        fields = []
        key_index = 0
        proxy_index = 0
        # id_field = -1
        coord_index = len(key_fields)
        velocity_index = len(key_fields)

        for key in key_fields:
            if coords_group and key in ['x', 'y', 'z']:
                if coords_group == 1:
                    coords_group = 2
                    coord_index = key_index
                    fields.append(numpy.zeros([frame_atoms, 3], dtype=float))
                    frame_info['format'].append('coord')
                    key_index += 1

                if key == 'x':
                    key_proxy[proxy_index] = [coord_index, 0]
                elif key == 'y':
                    key_proxy[proxy_index] = [coord_index, 1]
                elif key == 'z':
                    key_proxy[proxy_index] = [coord_index, 2]
                else:
                    log.error("Something strange happend"
                              " while assigning coordinates to field")

            elif coords_s_group and key in ['xs', 'ys', 'zs']:
                if coords_s_group == 1:
                    coords_s_group = 2
                    coord_index = key_index
                    fields.append(numpy.zeros([frame_atoms, 3], dtype=float))
                    # coord_scaled
                    frame_info['format'].append('coord_s')
                    key_index += 1

                if key == 'xs':
                    key_proxy[proxy_index] = [coord_index, 0]
                elif key == 'ys':
                    key_proxy[proxy_index] = [coord_index, 1]
                elif key == 'zs':
                    key_proxy[proxy_index] = [coord_index, 2]
                else:
                    log.error("Something strange happend"
                              " while assigning scaled coordinates to field")

            elif coords_u_group and key in ['xs', 'ys', 'zs']:
                if coords_u_group == 1:
                    coords_u_group = 2
                    coord_index = key_index
                    fields.append(numpy.zeros([frame_atoms, 3], dtype=float))
                    # coords unwrapped
                    frame_info['format'].append('coord_uw')
                    key_index += 1

                if key == 'xu':
                    key_proxy[proxy_index] = [coord_index, 0]
                elif key == 'yu':
                    key_proxy[proxy_index] = [coord_index, 1]
                elif key == 'zu':
                    key_proxy[proxy_index] = [coord_index, 2]
                else:
                    log.error("Something strange happend "
                              "while assigning unwrapped coordinates to field")

            elif velocity_group and key in ['vx', 'vy', 'vz']:
                if velocity_group == 1:
                    velocity_group = 2
                    velocity_index = key_index
                    fields.append(numpy.zeros([frame_atoms, 3], dtype=float))
                    frame_info['format'].append('velocity')
                    key_index += 1

                if key == 'vx':
                    key_proxy[proxy_index] = [velocity_index, 0]
                elif key == 'vy':
                    key_proxy[proxy_index] = [velocity_index, 1]
                elif key == 'vz':
                    key_proxy[proxy_index] = [velocity_index, 2]
                else:
                    log.error("Something strange happend"
                              " while assigning fields")

            elif key == 'type':
                fields.append(numpy.zeros(frame_atoms, dtype=int))
                frame_info['format'].append(key)
                key_proxy[proxy_index] = key_index
                key_index += 1

            elif key == 'id':
                # id_field = key_index
                fields.append(numpy.zeros(frame_atoms, dtype=int))
                frame_info['format'].append(key)
                key_proxy[proxy_index] = key_index
                key_index += 1

            elif key == 'element':
                frame_info['format'].append(key)
                # *TO DO* upgrade from simple python string list
                fields.append([None]*frame_atoms)
                key_proxy[proxy_index] = key_index
                key_index += 1
            else:
                fields.append(numpy.zeros(frame_atoms, dtype=float))
                frame_info['format'].append(key)
                key_proxy[proxy_index] = key_index
                key_index += 1
            proxy_index += 1

        if skip:
            for i in xrange(frame_atoms):
                infile.readline()
            return 1

        fields_n = proxy_index
        for i in xrange(frame_atoms):
            sp = infile.readline().split()

            for j in xrange(fields_n):
                jproxy = key_proxy[j]

                if isinstance(jproxy, list):
                    fields[jproxy[0]][i][jproxy[1]] = sp[j]
                elif jproxy is not None:
                    fields[jproxy][i] = sp[j]

        frame = tools.Frame()

        for key in frame_info.keys():
            frame.add(key, frame_info[key])

        for i in range(len(fields)):
            frame.add(frame_info['format'][i], fields[i])


        return frame

    def _read_xyz(self, skip):
        """
        ---
        """
        frame = tools.Frame()
        infile = self.infile

        line = infile.readline()
        if not line:
            return False

        frame_atoms = int(line)
        frame.add('size', int(frame_atoms))

        line = infile.readline()
        frame.add('comment', line)
        frame.add('format', ['element', 'coord'])

        coords = numpy.empty([frame_atoms, 3], dtype=float)
        elements = [None]*frame_atoms

        if skip:
            for i in xrange(frame_atoms):
                infile.readline()
            return 1

        for i in xrange(frame_atoms):
            sl = infile.readline().split()
            elements[i] = sl[0]
            coords[i][0] = sl[1]
            coords[i][1] = sl[2]
            coords[i][2] = sl[3]
        frame.add('element', elements)
        frame.add('coord', coords)
        return frame


class RMS:
    """
    Tool for calculating roughness of atomic scale surfaces

    Args:
        coords (numpy float list (n, 3)):
        r_sample : vdw radius of sample atoms
        r_probe : vdw radius of probe atom
        sampling_box (optional array (2,2)): box in which sampling array/matrix
        will be created
        active_box (optional array (3,2)): box in which the atoms are included
        in rms calcuation, if None - program will find optimal box
        offset : for sampling array probing (final offset = offset * r_vdw)
    """
    def __init__(self, coords, r_sample, r_probe, sampling_box=None,
                 active_box=None, offset=1.0):
        self.coords = coords         # array of surface atoms coordinates
        self.r_sample = r_sample     # r vdw of surface atoms
        self.r_probe = r_probe       # r vdw of probe atom
        self.r = r_sample + r_probe  # total bond r vdw surface-probe
        self.r2 = self.r*self.r
        self.offset = offset*self.r

        self.sampling_array_flag = False
        self.sampling_array_dim_flag = False
        self.active_box_flag = False
        self.surf_array_flag = False

        # set or find box in which surface atoms will be searched
        if active_box:
            self.set_active_box(active_box)
        else:
            self.find_activ_box()

        if sampling_box:
            self.set_sampling_box(sampling_box)

    def set_active_box(self, active_box):
        """
        Set up active_box variable after checking correctness of the argument

        Args:
            active_box (list): [[xmin, xmax], [ymin, ymax], [zmin, zmax]]
                only atoms in box will be taken into account during analysis
        """
        try:
            self.active_box = numpy.array([[active_box[0][0], active_box[0][1]],
                                           [active_box[1][0], active_box[1][1]],
                                           [active_box[2][0], active_box[2][1]]],
                                          dtype=float)
            log.info("Active box set:\n" + str(self.active_box))

        except:
            log.error("Wrong active_box variable format, should be:"
                      " [[xmin, xmax], [ymin, ymax], [zmin, zmax]]")

        if (active_box[0][0] >= active_box[0][1]):
            log.error("Ill-defined X demension for box")
        if (active_box[1][0] >= active_box[1][1]):
            log.error("Ill-defined Y demension for box")
        if (active_box[2][0] >= active_box[2][1]):
            log.error("Ill-defined Z demension for box")

        self.active_box_flag = True

    def find_activ_box(self, h_ratio=1):
        """
        Find box large enough to fit all atoms

        Args:
            h_ratio : what part of the atoms in the box will be used for
            claculations (atoms with higher z component will be accepted)
        """

        coords = self.coords
        xmax, ymax, zmax = coords.max(axis=0)
        xmin, ymin, zmin = coords.min(axis=0)

        self.set_active_box([[xmin, xmax], [ymin, ymax], [zmin, zmax]])

    def set_sampling_box(self, sampling_box):
        """
        Set up active_box variable after checking correctness of the argument

        Args:
            active_box (list): [[xmin, xmax], [ymin, ymax], [zmin, zmax]]
                only atoms in box will be taken into account during analysis
        """
        try:
            self.sampling_box = numpy.array([[sampling_box[0][0], sampling_box[0][1]],
                                             [sampling_box[1][0], sampling_box[1][1]]],
                                             dtype=float)
        except:
            log.error("Wrong sampling_box variable format, should be:"
                      " [[xmin, xmax], [ymin, ymax]]")

        if (sampling_box[0][0] >= sampling_box[0][1]):
            log.error("Ill-defined X demension for box")
        if (sampling_box[1][0] >= sampling_box[1][1]):
            log.error("Ill-defined Y demension for box")

    def get_surf_array_oa2d(self, lc_rep, format='array'):
        """
        oa2d - one atom (per cell) 2 dimensional array

        Assigns atoms to 2D grid in xy plane. Each created cell contains
        information only about atom with the largest z.

        Larger lc_rep results in improved accuracy at the cost of computing time

        Caution:
            array may contain nan value

        Args:
            sampling_array (numpy.ndarray): array of point for which height
                will be calculated
            lc_rep (int): how many logic cells will be taken to calculate
                one point 1 - (square 3x3) 2 - (square 5x5) 3 - (square 7x7)...
            format (str): return array format
                array - (3, N) coords array
                matrix - (N, N) matrix of z value
        """
        if format not in ['array', 'matrix']:
            log.error("Unknow surf array format: %s" % str(format))

        if not self.active_box_flag:
            self.find_activ_box()

        if not self.sampling_array_flag:
            self.calc_sampling_array()

        if format == 'matrix' and not self.sampling_array_dim_flag:
            log.error("format=matrix can only be used with automatically"
                      "build surface_array")
            return None

        # Step 1 - create array of logica cells
        # ------------------------------------
        lc_r = self.r/float(lc_rep)
        lc_rep = int(lc_rep)
        active_box = self.active_box
        sampling_array = self.sampling_array

        # + 0.01 is there to make sure the boundary conditions are no violated
        x_dimension = active_box[0][1] - active_box[0][0] + 0.01
        y_dimension = active_box[1][1] - active_box[1][0] + 0.01

        x_rep = int(math.ceil(x_dimension/lc_r))
        y_rep = int(math.ceil(y_dimension/lc_r))

        logical_cells = numpy.zeros((x_rep, y_rep, 3), dtype=float)
        logical_cells.fill(active_box[2][0])

        x_shift = -active_box[0][0]
        y_shift = -active_box[1][0]

        coords = self.coords
        # loop over all atoms and assign them to the appropriate cells
        for i in xrange(len(coords)):
            x_index = int(math.floor(coords[i][0]+x_shift)/lc_r)
            y_index = int(math.floor(coords[i][1]+y_shift)/lc_r)
            z = coords[i][2]

            # only one atom with highest z per logic cell
            if z > logical_cells[x_index][y_index][2]:
                logical_cells[x_index][y_index][0] = coords[i][0]
                logical_cells[x_index][y_index][1] = coords[i][1]
                logical_cells[x_index][y_index][2] = z

        # Step 2 - using generated logical cells and points to sample
        # create array with corresponding heights
        # ------------------------------------

        x_index_max = int(x_rep - lc_rep)
        y_index_max = int(y_rep - lc_rep)

        # neigh_template added to cell position gives neighbors of that cell
        neigh_template = numpy.arange(-lc_rep, lc_rep+1, 1, dtype=int)

        surface_array = numpy.empty((len(sampling_array), 3), dtype=float)
        for i in xrange(len(sampling_array)):
            # calculate index  corresponding to probe x y
            probe_xy = sampling_array[i]
            x_index = int(math.floor(probe_xy[0]+x_shift)/lc_r)
            y_index = int(math.floor(probe_xy[1]+y_shift)/lc_r)

            # check boundary conditions
            if x_index < lc_rep:
                x_neigh = neigh_template[(lc_rep-x_index):]+x_index
            elif x_index >= x_index_max:
                x_neigh = neigh_template[:(lc_rep+x_rep-x_index)]+x_index
            else:
                x_neigh = neigh_template+x_index

            if y_index < lc_rep:
                y_neigh = neigh_template[(lc_rep-y_index):]+y_index
            elif y_index >= y_index_max:
                y_neigh = neigh_template[:(lc_rep+y_rep-y_index)]+y_index
            else:
                y_neigh = neigh_template+y_index

            z = None
            # loop over all neighbors
            for x_n_index in x_neigh:
                for y_n_index in y_neigh:
                    z_new = self.calc_z_oa2d(probe_xy,
                                         logical_cells[x_n_index][y_n_index])
                    if z_new > z:
                        z = z_new
            surface_array[i][0] = probe_xy[0]
            surface_array[i][1] = probe_xy[1]
            surface_array[i][2] = z

        if format == 'array':
            self.surf_array_flag = True
            return surface_array

        # map array on matrix
        if format == 'matrix':
            xdim = self.sampling_array_x
            ydim = self.sampling_array_y
            surface_matrix = numpy.empty((xdim, ydim), dtype=float)
            for xi in xrange(xdim):
                for yi in xrange(ydim):
                    totali = xi*ydim + yi
                    surface_matrix[xi][yi] = surface_array[totali][2]
            return surface_matrix

    def calc_z_oa2d(self, probe_xy, sample_atom):
        """
        For a give pair of atoms probe-sample, calculate the maximum z which
        probe atom can achive.

        Args:
            probe_xy (numpy.ndarray):
            sample_atom (numpy.ndarray):

        Returns:
            None: if probe atom and sample atom does not touch
            float: maximu z of probe atom
        """
        r2 = self.r2
        dx = probe_xy[0] - sample_atom[0]
        dy = probe_xy[1] - sample_atom[1]
        dxy2 = dx*dx + dy*dy
        z = sample_atom[2]

        if dxy2 > r2:
            return None
        else:
            return math.sqrt(r2-dxy2)+z

    def calc_sampling_array(self, jump=0.5):
        """
        Generate and set evenly spaced sampling array in xy plane

        Args:
            jump (optional float): grid spacing
        """
        offset = self.offset
        if not self.active_box_flag:
            self.find_activ_box()

        dim_x = self.sampling_box[0]
        dim_y = self.sampling_box[1]

        x = numpy.arange(dim_x[0] + offset, dim_x[1] - offset, jump)
        y = numpy.arange(dim_y[0] + offset, dim_y[1] - offset, jump)

        self.sampling_array_x = len(x)
        self.sampling_array_y = len(y)
        self.sampling_array_dim_flag = True

        X, Y = numpy.meshgrid(x, y)
        self.sampling_array = numpy.array([Y.flatten(), X.flatten()]).T
        self.sampling_array_flag = True

    def set_sampling_array(self, array):
        """
        Sets sampling poins

        Args:
            array (numpy.ndarray): arry in xy plane of sampling points
        """

        self.sampling_array = array
        self.sampling_array_flag = True

    def compute(self, lc_rep=2):
        """
        Compute rms for previously passed parameters

        Args:
            lc_rep (int): how many logic cells will be taken to calculate one
            point, 2 - (square 3x3) 2 - (square 5x5) 3 - (square 7x7)  ...

        float: rms
        """
        if not self.surf_array_flag:
            surf_array = self.get_surf_array_oa2d(lc_rep, format='array')

        z_ave = numpy.average(surf_array[:,2])
        z_sq_sum = 0.0

        valid_cnt = 0
        for i in xrange(len(surf_array)):
            z_current = surf_array[i][2]
            if numpy.isnan(z_current):
                continue
            valid_cnt += 1
            res = z_ave-surf_array[i][2]
            z_sq_sum += res*res

        return math.sqrt(z_sq_sum/valid_cnt)

