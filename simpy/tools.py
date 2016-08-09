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
