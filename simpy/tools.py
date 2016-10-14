import numpy as np
import logging
log = logging.getLogger(__name__)

class BinAtom:
    def __init__(self, size, bin_size):
        self.size = size
        self.bin_size = bin_size

        self.init_bins()

    def init_bins(self):
        bins_dim = []

        for dim in self.size:
            bins_num = int(np.ceil(float(dim[1] - dim[0])/self.bin_size))
            bins_dim.append(bins_num)

        self.bins_dim = bins_dim
        self.bins = np.zeros((bins_dim[0], bins_dim[1], bins_dim[2]),
                             dtype=float)

    def add(self, frame, weight=None):
        size = self.size

        coord = frame['coord']
        selx = np.logical_and(coord[:,0] > size[0][0], coord[:,0] < size[0][1])
        sely = np.logical_and(coord[:,1] > size[1][0], coord[:,1] < size[1][1])
        selz = np.logical_and(coord[:,2] > size[2][0], coord[:,2] < size[2][1])

        sel = np.logical_and(np.logical_and(selz, sely), selx)

        log.info("Selected: %i / %i" % (np.sum(sel), len(sel)))

        # slow implementation (without digitize form np)
        selected = coord[sel]

        size = self.size
        bin_size = self.bin_size
        bins = self.bins
        start_pos = [size[0][0], size[1][0], size[2][0]]
        print start_pos
        for i in xrange(len(selected)):
            acoord = selected[i]

            xi = int((acoord[0] - start_pos[0])/bin_size)
            yi = int((acoord[1] - start_pos[1])/bin_size)
            zi = int((acoord[2] - start_pos[2])/bin_size)

            bins[xi][yi][zi] += 1.0;


    def project_bins(self, axis=2, center=0.0, wmethod="linear"):
        print "reutrn"
        return np.sum(self.bins, axis=axis)


    def as_frame(self, accept=1.0):
        bins = self.bins
        bins_dim = self.bins_dim

        for xi in  xrange(bins_dim[0]):
            for yi in  xrange(bins_dim[1]):
                for zi in  xrange(bins_dim[2]):
                    if(bins[xi][yi][zi] >= accept):
                        pass




class Frame:
    def __init__(self):
        self.record = {}

    def __getitem__(self, key):
        if key in self.record:
            return self.record[key]
        return None

    def __len__(self):
        if 'coord' in self.record:
            return len(self.record['coord'])
        return None


    def keys(self):
        return self.record.keys()

    def add(self, key, value):
        self.record[key] = value

    def gen_types(self, atom_type=1):
        atoms = len(self['coord'])
        self.add('type', [atom_type]*atoms)

    def gen_element(self, element_name='H'):
        atoms = len(self['coord'])
        self.add('element', [element_name]*atoms)

    def select(self, box):
        """
        return frame containing only records inside box
        box=[[xmin, xmax], [ymin, ymax], [zmin, zmax]]

        current:
            only coordinates
        """

        coord = self['coord']
        minCr = np.array([box[0][0], box[1][0], box[2][0]])
        maxCr = np.array([box[0][1], box[1][1], box[2][1]])

        inside = np.all(np.logical_and(minCr <= coord, coord <= maxCr), axis=1)
        newCoord = coord[inside]

        oFrame = Frame()
        oFrame.add('coord', newCoord)
        return oFrame
