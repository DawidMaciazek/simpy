"""
This module is used to analyze simulation results.



"""
import numpy
# TO DO: set proper logger! !!
import logging as log
log.basicConfig(format="%(levelname)s: %(message)s", level=10)


class traj:
    def __init__(self, filename, format='lammpstrj', log_level=10):
        """
        Args:
            filename: filename for reading
            format: format of the readed file
        """

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

    def read(self, read_n=-1):
        """
        Reads trajectory file frame by frame

        Args:
            read_n: number of trajestories to read
            read_n = -1: read the whole file
        Returns:
            list of frames
        """

        read_whole = False if read_n > 0 else True
        frame = True

        frames = []
        log.info("Read %i frames" % read_n)
        frames_readed = 0
        read_begining = str(read_n) if read_n > 0 else "INF"
        while frame and (read_n or read_whole):
            frame = self.call_read(False)
            if frame:
                frames_readed += 1
                frame.append(frame)
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
                    frame_info['format'].append('coords')
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
                    frame_info['format'].append('coords_scaled')
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
                    frame_info['format'].append('coords_unwrapped')
                    key_index += 1

                if key == 'xu':
                    key_proxy[proxy_index] = [coord_index, 0]
                elif key == 'yu':
                    key_proxy[proxy_index] = [coord_index, 1]
                elif key == 'zu':
                    key_proxy[proxy_index] = [coord_index, 2]
                else:
                    log.error("Something strange happend"
                              " while assigning unwrapped coordinates to field")

            elif velocity_group and key in ['vx', 'vy', 'vz']:
                if velocity_group == 1:
                    velocity_group = 2
                    velocity_index = key_index
                    fields.append(numpy.zeros([frame_atoms, 3], dtype=float))
                    frame_info['velocity'].append('velocity')
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

        frame_output = [frame_info]
        frame_output.extend(fields)
        return frame_output

    def _read_xyz(self, skip):

        infile = self.infile
        frame_info = {}
        line = infile.readline()
        if not line:
            return False

        frame_atoms = int(line)
        frame_info['atoms'] = int(frame_atoms)

        line = infile.readline()
        frame_info['comment'] = line
        # * TO DO * allow the acceptance of dditional fields
        frame_info['fields'] = ['element', 'coords']

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
        return [frame_info, elements, coords]

class rms:
    def __init__(self, coords, active_box=None):
        # basic plan:
        # 1 chek avai area
        # 2 create 2D height map
        # 3 calculate rms
        #



        if active_box:
            self.set_active_box(active_box)
            pass
        else:
            self.find_activ_box()


    def set_active_box(self, active_box):
        try:
            self.active_box = numpy([active_box[0][0], active_box[0][1],
                                    active_box[1][0], active_box[1][1],
                                    active_box[2][0], active_box[2][1]])
        except:
            log.error("Wrong active_box variable format, should be:"
                        " [xmin, xmax], [ymin, ymax], [zmin, zmax]")


        # * TO DA * check if *min < *max

    def find_activ_box(self, h_ratio=0.5):
        coords = self.coords
        xmax, ymax, zmax = coords.max(axis=0)
        xmin, ymin, zmin = coords.min(axis=0)



