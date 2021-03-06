import numpy as np
import logging
import copy
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
        for i in xrange(len(selected)):
            acoord = selected[i]

            xi = int((acoord[0] - start_pos[0])/bin_size)
            yi = int((acoord[1] - start_pos[1])/bin_size)
            zi = int((acoord[2] - start_pos[2])/bin_size)

            bins[xi][yi][zi] += 1.0;


    def project_bins(self, axis=2, center=0.0, wmethod="linear"):
        return np.sum(self.bins, axis=axis)


    def as_frame(self, accept=1.0):
        bins = self.bins
        bins_dim = self.bins_dim

        for xi in  xrange(bins_dim[0]):
            for yi in  xrange(bins_dim[1]):
                for zi in  xrange(bins_dim[2]):
                    if(bins[xi][yi][zi] >= accept):
                        pass


def digitize3d(frame, boundary, bin_size):
    #frame = frame.select(boundary)
    boundary = np.array(boundary, dtype=float)
    bin_steps = (boundary[:,1]-boundary[:,0])/bin_size + 1
    bin_centers = [None] * 3
    bin_edges = [None] * 3
    bin_digitized = [None] * 3
    for i in range(3):
        bin_centers[i] = np.linspace(boundary[i][0], boundary[i][1], bin_steps[i])
        bin_edges[i] = bin_centers[i][:-1]+(0.5*bin_size)
        bin_digitized[i] = np.digitize(frame['coord'][:,i], bin_edges[i])
    #bin_centers = np.array(bin_centers, dtype=float)
    #bin_edges = np.array(bin_edges, dtype=float)
    bin_digitized = np.array(bin_digitized, dtype=int).T

    return  bin_digitized, bin_centers



class Frame:
    def __init__(self, coord=None, element=None):
        self.record = {}
        if coord is not None:
            self.add("coord", coord)
        if element is not None:
            self.add("element", element)


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
        atoms = self.__len__()
        self.add('type', [atom_type]*atoms)

    def gen_element(self, element_name='H'):
        atoms = self.__len__()
        self.add('element', [element_name]*atoms)

    def types_to_element(self, edict):
        self.add('element', np.vectorize(edict.get)(self.record['type']))

    def remove_bool_list(self, remove_list):
        keep_list = np.logical_not(remove_list)
        removed_atoms = sum(remove_list)
        for key in self['format']:
            self.record[key] = self[key][keep_list]
            log.info("{} atoms removed from field {}".format(removed_atoms, key))
        # update len
        self.record["atoms"] = self.__len__()

    def remove_element(self, element_name):
        if "element" not in self.keys():
            log.warning("This frame does not have element field - nothing will be removed")
            return

        log.info("Removing {} element".format(element_name))

        elements = self['element']
        remove_list = (elements == element_name)
        self.remove_bool_list(remove_list)

    def remove_type(self, type_name):
        if "type" not in self.keys():
            log.warning("This frame does not have type field - nothing will be removed")
            return

        log.info("Removing {} type".format(type_name))

        type = self['type']
        remove_list = (type == type_name)
        self.remove_bool_list(remove_list)

    def remove_lesser_than(self, value, axis):
        if 'coord' not in self.record:
            log.warn("This frame does not have suitable coordinate system")
            return

        coords = self.record['coord']
        remove_list = coords[:,axis] < value
        self.remove_bool_list(remove_list)

    def periodic(self, period):
        if 'box' not in self.keys():
            log.warning("Periodic operation requires box field")
            return
        box = self.record['box']

        box_size = []
        box_size.append(box[0][1] - box[0][0])
        box_size.append(box[1][1] - box[1][0])
        box_size.append(box[2][1] - box[2][0])

        for i in range(3):
            base = {}
            for key in self.keys():
                base[key] = self.record[key]

            for j in range(period[i][0], period[i][1]+1):
                if j == 0:
                    continue

                shifted_base = np.zeros(base["coord"].shape, dtype=float)
                shifted_base[:,i] += j*box_size[i]

                self.record["coord"] = np.append(self.record["coord"], base["coord"] + shifted_base, axis=0)

                for key in self.record["format"]:
                    if key == "coord":
                        continue
                    self.record[key] = np.append(self.record[key], base[key])

        # update box size and atoms
        self.record["atoms"] = self.__len__()
        self.record["box"] = [ [box[i][0] + box_size[i]*period[i][0], box[i][1] + box_size[i]*period[i][1]] for i in range(3) ]


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
        try:
            oFrame.add("element", self["element"])
        except:
            pass
        return oFrame

    def copy(self):
        cframe = Frame()
        for key in self.record:
            try:
                if isinstance(self.record[key], list):
                    cframe.add(key, copy.deepcopy(self.record[key]))
                else:
                    cframe.add(key, self.record[key].copy())
            except AttributeError:
                cframe.add(key, self.record[key])

        return cframe
