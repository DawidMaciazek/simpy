class Frame:
    def __init__(self):
        self.record = {}

    def add(self, key, value):
        self.record[key] = value

    def __getitem__(self, key):
        return self.record[key]
