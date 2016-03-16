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
        while frame and (read_n or read_whole):
            print "current %i" % read_n
            frame = self.call_read(False)
            frames.append(frame)
            read_n -= 1

        return frames

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
                timestep = int(infile.readline()[0])
                frame_info['timestep'] = int(timestep)

            elif key_field_name == "NUMBER":
                frame_atoms = int(infile.readline()[0])
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

        if 'vx' in key_fields and 'vy' in key_fields and 'vz' in key_fields:
            velocity_group = 1
        else:
            velocity_group = False

        addInfo['format'] = []

        # transform key position index in file to index in output list
        key_proxy = [None]*len(key_fields)

        fields = []
        key_index = 0
        proxy_index = 0
        id_field = -1
        coord_index = len(key_fields)
        velocity_index = len(key_fields)

        for key in key_fields:
            if coords_group and key in ['x', 'y', 'z']:
                if coords_group == 1:
                    coord_index = key_index
                    coords_group = 2
                    fields.append(numpy.zeros([frame_atoms, 3], dtype=float))
                if key == 'x':
                    key_proxy[proxy_index] = [coord_index, 0]
                elif key == 'y':
                    key_proxy[proxy_index] = [coord_index, 1]
                elif key == 'z':
                    key_proxy[proxy_index] = [coord_index, 2]
                else:
                    log.error("Something strange happend while assigning fields")

            if velocity_group and key in ['vx', 'vy', 'vz']:
                if velocity_group == 1:
                    velocity_group = 2
                    fields.append(numpy.zeros([frame_atoms, 3], dtype=float))

                if key == 'vx':
                    key_proxy[proxy_index] = [velocity_index, 0]
                elif key == 'vy':
                    key_proxy[proxy_index] = [velocity_index, 1]
                elif key == 'vz':
                    key_proxy[proxy_index] = [velocity_index, 2]
                else:
                    log.error("Something strange happend while assigning fields")

            elif key == 'type':
                pass
            elif key == 'id':
                pass
            elif key == 'element':
                pass
            else:
                pass



    def _read_xyz(self, skip):
        print "calling XYZ ! function)"
