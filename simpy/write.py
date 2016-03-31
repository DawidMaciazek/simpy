class xyz:
    def __init__(self, filename):
        self.outfile = open(filename, 'w')

    def write(self, types, atoms, optiona=None, comment="generated by sympy lib"):
        outfile = self.outfile
        atoms_n = len(atoms)

        outfile.write("%i\n%s\n" % (atoms_n, str(comment)))
        for i in xrange(atoms_n):
            outfile.write("%s %f %f %f\n" %
                          (types[i], atoms[i][0], atoms[i][1], atoms[i][2]))

    def close(self):
        self.outfile.close()
