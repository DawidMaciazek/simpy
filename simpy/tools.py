import numpy as np

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
