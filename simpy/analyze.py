"""
This module is used to analyze simulation results.



"""

class traj:
    def __init__(self, filename, format='lammpstrj'):
        """
        Args:
            filename: filename for reading
            format:
        """

        function_name = "_read_" + format

        try:
            self.call_read = getattr(self, function_name)
        except AttributeError:
            print "Unknown format: %s" % str(format)


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
        print "calling lammpstrj function)"

        return True

    def _read_xyz(self, skip):
        print "calling XYZ ! function)"
