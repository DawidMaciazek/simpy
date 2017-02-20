import numpy
import math
import math_tools

import logging
log = logging.getLogger(__name__)

def check_if_numpy2D(array, fun_name):
    # check if numpy array
    if not type(array).__module__ == numpy.__name__:
        log.error("function [%s] only accepts numpy array" % fun_name)

    if len(array.shape) != 2 or array.shape[1] != 3:
        log.error("function [%s] only accepts numpy array of shape (n, 3)"
                  % fun_name)

class Lattice:
    def __init__(self, box_size=None, lattice_system=None, origin=None,
                 lat_constant=None, miller=None, miller_angle=None,
                 lat_vectors=None, ):
        # Setting up default parameters
        if box_size:
            self.set_box_size(box_size)
        else:
            self.set_box_size(
                numpy.array([[-10., 10.], [-10., 10.], [-10., 10.]]))

        if lattice_system:
            self.set_system(lattice_system)
        else:
            self.set_system("pcc")

        if origin:
            self.set_origin(origin)
        else:
            self.set_origin(numpy.array([0., 0., 0.]))

        if lat_vectors:
            self.set_lattice_vectors(lat_vectors)
        else:
            self.set_lattice_vectors(numpy.array([[2., 0., 0.], [0., 2., 0.], [0., 0., 2.]]))

        if lat_constant:
            if not miller_angle:
                miller_angle=0.0

            if not miller:
                miller=[1,0,0]

            self.set_lattice_miller(lat_constant, miller, miller_angle)

    def set_box_size(self, box_size):
        try:
            xp = numpy.array([float(xp_) for xp_ in box_size[0]])
            yp = numpy.array([float(yp_) for yp_ in box_size[1]])
            zp = numpy.array([float(zp_) for zp_ in box_size[2]])
        except:
            log.error("Wrong box size variable format, should be:\n"
                      "[[x_min, x_max], [ymin, ymax], [zmin, zmax]")

        if (xp[0] >= xp[1]):
            log.error("Ill-defined X demension for box")
        if (yp[0] >= yp[1]):
            log.error("Ill-defined Y demension for box")
        if (zp[0] >= zp[1]):
            log.error("Ill-defined Z demension for box")

        self.box_size = numpy.array([xp, yp, zp], dtype=float)
        log.info("Current box size:\n%s" % str(self.box_size))

        # set up box corners
        bs = self.box_size
        self.box_corners = numpy.array([
            [bs[0][0], bs[1][0], bs[2][0]],  # 0
            [bs[0][0], bs[1][0], bs[2][1]],  # 1
            [bs[0][0], bs[1][1], bs[2][0]],  # 2
            [bs[0][0], bs[1][1], bs[2][1]],  # 3
            [bs[0][1], bs[1][0], bs[2][0]],  # 4
            [bs[0][1], bs[1][0], bs[2][1]],  # 5
            [bs[0][1], bs[1][1], bs[2][0]],  # 6
            [bs[0][1], bs[1][1], bs[2][1]],  # 7
        ])

    def set_system(self, lat_type):
        # simple
        if lat_type in ["pcc", "simple", "P", "p"]:
            self.lattice_system = "pcc"

        # Body Centered Cubic
        elif lat_type in ["bcc", "B", "b"]:
            self.lattice_system = "bcc"

        # Face Centered Cubic
        elif lat_type in ["fcc", "F", "f"]:
            self.lattice_system = "fcc"

        # Diamond
        elif lat_type in ["diamond", "D", "d"]:
            self.lattice_system = "diamond"

        else:
            log.error("Unknow lattice type: %s\n"
                      "available types: pcc, bcc, fcc, (?daimond?)" % lat_type)

        log.info("Current lattice system: %s" % lat_type)

    def set_origin(self, origin):
        try:
            self.origin = numpy.array([float(origin[0]), float(origin[1]),
                                      float(origin[2])], dtype=float)
        except:
            log.error("Wrong origin variable format, should be:\n"
                      "[x_origin, y_origin, z_origin]")

        log.info("Current origin: %s" % str(self.origin))

    def set_lattice_vectors(self, lat_vectors):
        try:
            al = numpy.array([float(ap_) for ap_ in lat_vectors[0]], dtype=float)
            bl = numpy.array([float(bp_) for bp_ in lat_vectors[1]], dtype=float)
            cl = numpy.array([float(cp_) for cp_ in lat_vectors[2]], dtype=float)
        except:
            log.error("Wrong lat_vectors format, should be:\n"
                      "[[ax, ay, az], [bx, by, bz], [cx, cy, cz]]")

        if(numpy.dot(al, al) == 0):
            log.error("A vector has zero lenght")
        if(numpy.dot(bl, bl) == 0):
            log.error("B vector has zero lenght")
        if(numpy.dot(cl, cl) == 0):
            log.error("C vector has zero lenght")

        cross_AB = numpy.cross(al, bl)
        if(numpy.dot(cross_AB, cross_AB) == 0):
            log.error("A x B = 0, vectors may by parallel")

        if(numpy.dot(cross_AB, cl) == 0):
            log.error("C lies on AB plane")

        if(numpy.dot(cross_AB, cl) < 0):
            log.error("A B C vectros don't form right-handed system")

        self.lattice_vectors = numpy.array([al, bl, cl], dtype=float)
        log.info("Current lattice vectors:\n%s" % str(self.lattice_vectors))

    def set_lattice_miller(self, lat_constant, miller, miller_angle):
        log.info("Using lattice constant: {}, miller vector: {} and miller angle {}".format(lat_constant, miller,miller_angle))

        try:
            miller = numpy.array([miller[0], miller[1], miller[2]], dtype=float)
        except:
            log.error("Bad Miller indexes format: %s\n shoud be [a, b, c]" %
                      str(miller))

        self.set_lattice_vectors([[lat_constant, 0., 0.],
                                  [0., lat_constant, 0.],
                                  [0., 0., lat_constant]])
        self.set_miller(miller, miller_angle)

    def set_miller(self, miller, miller_angle=0.0, orient=[0., 0., 1.]):
        lattice_vectors = self.lattice_vectors
        direction = numpy.array(orient, dtype=float)
        miller = numpy.array(miller)

        # normal vector to the desired miller plane (miller vector)
        lattice_vectors_sum = numpy.sum(lattice_vectors, axis=0)
        miller = numpy.array([miller[0], miller[1],
                              miller[2]]) * lattice_vectors_sum

        rotation_vector = numpy.cross(miller, direction)
        rotation_vector = rotation_vector / numpy.linalg.norm(rotation_vector)

        rotation_angle = math.acos(numpy.dot(miller, direction) /
                                   (numpy.linalg.norm(miller) *
                                    numpy.linalg.norm(direction)))

        rmatrix = math_tools.rotation_matrix(rotation_vector, rotation_angle)
        rmatrix_main = math_tools.rotation_matrix(direction, miller_angle)

        lattice_final = [None]*3
        for i in xrange(3):
            lattice_final[i] = numpy.dot(rmatrix_main,
                                         numpy.dot(rmatrix, lattice_vectors[i]))

        self.set_lattice_vectors(numpy.array(lattice_final))

    def get_nodes(self):
        a = self.lattice_vectors[0]
        b = self.lattice_vectors[1]
        c = self.lattice_vectors[2]

        # calculate jumps size and starting point
        plane_point = numpy.empty((3, 3), dtype=float)
        norm_vec = numpy.empty((3, 3), dtype=float)

        # number of jumps in each dimension
        jumps = numpy.empty(3, dtype=int)

        jumps[0], norm_vec[0], plane_point[0] = self.find_jumps_vec(b, c, a)
        jumps[1], norm_vec[1], plane_point[1] = self.find_jumps_vec(c, a, b)
        jumps[2], norm_vec[2], plane_point[2] = self.find_jumps_vec(a, b, c)

        plane_factor = numpy.empty(3, dtype=float)

        for i in xrange(3):
            plane_factor[i] = (norm_vec[i][0]*plane_point[i][0] +
                               norm_vec[i][1]*plane_point[i][1] +
                               norm_vec[i][2]*plane_point[i][2])

        # find planes cross
        start_vec = numpy.linalg.solve(norm_vec, plane_factor)

        size_max = (int(jumps[0])+1)*(int(jumps[1])+1)*(int(jumps[2])+1)
        nodes = numpy.empty((size_max, 3), dtype=float)
        node_id = 0

        clvl = numpy.copy(start_vec)
        for ci in xrange(int(jumps[2])+1):
            blvl = numpy.copy(clvl)
            for bi in xrange(int(jumps[1])+1):
                alvl = numpy.copy(blvl)
                for ai in xrange(int(jumps[0])+1):
                    nodes[node_id][0] = alvl[0]
                    nodes[node_id][1] = alvl[1]
                    nodes[node_id][2] = alvl[2]
                    node_id += 1
                    alvl += a
                blvl += b
            clvl += c

        lattice_vectors = self.lattice_vectors
        bs = self.box_size
        if (self.lattice_system == "pcc"):
            pass

        elif (self.lattice_system == "bcc"):
            nodes_ex = nodes + numpy.sum(lattice_vectors * 0.5, axis=0)
            nodes = numpy.concatenate((nodes, nodes_ex), axis=0)

        elif (self.lattice_system == "fcc"):
            nodes_ex0 = nodes + (lattice_vectors[0] + lattice_vectors[1]) * 0.5
            nodes_ex1 = nodes + (lattice_vectors[1] + lattice_vectors[2]) * 0.5
            nodes_ex2 = nodes + (lattice_vectors[2] + lattice_vectors[0]) * 0.5

            nodes = numpy.concatenate((nodes, nodes_ex0, nodes_ex1, nodes_ex2),
                                      axis=0)

        elif (self.lattice_system == "diamond"):
            nodes_ex0 = nodes + (lattice_vectors[0] + lattice_vectors[1]) * 0.5
            nodes_ex1 = nodes + (lattice_vectors[1] + lattice_vectors[2]) * 0.5
            nodes_ex2 = nodes + (lattice_vectors[2] + lattice_vectors[0]) * 0.5

            nodes_ex3 = nodes + numpy.sum(lattice_vectors * 0.25, axis=0)
            nodes_ex4 = nodes + (lattice_vectors[0]*0.75 +
                                 lattice_vectors[1]*0.75 +
                                 lattice_vectors[2]*0.25)
            nodes_ex5 = nodes + (lattice_vectors[0]*0.25 +
                                 lattice_vectors[1]*0.75 +
                                 lattice_vectors[2]*0.75)
            nodes_ex6 = nodes + (lattice_vectors[0]*0.75 +
                                 lattice_vectors[1]*0.25 +
                                 lattice_vectors[2]*0.75)

            nodes = numpy.concatenate((nodes, nodes_ex0, nodes_ex1, nodes_ex2,
                                       nodes_ex3, nodes_ex4, nodes_ex5,
                                       nodes_ex6), axis=0)

        x_list = numpy.logical_and(nodes[:, 0] > bs[0][0],
                                   nodes[:, 0] < bs[0][1])
        y_list = numpy.logical_and(nodes[:, 1] > bs[1][0],
                                   nodes[:, 1] < bs[1][1])
        z_list = numpy.logical_and(nodes[:, 2] > bs[2][0],
                                   nodes[:, 2] < bs[2][1])

        final_list = numpy.logical_and(z_list,
                                       numpy.logical_and(x_list, y_list))

        nodes = nodes[final_list]
        return nodes


    def find_jumps_vec(self, v1, v2, v3):
        # calculate normal vector to the plane
        norm_vec = numpy.cross(v1, v2)
        norm_vec = norm_vec/numpy.linalg.norm(norm_vec)

        # find max and min vectors
        corners = self.box_corners

        mincr = numpy.dot(corners[0], norm_vec)
        maxcr = numpy.dot(corners[0], norm_vec)
        for crvec in corners:
            crlen = numpy.dot(crvec, norm_vec)
            if(crlen < mincr):
                mincr = crlen
            elif(crlen > maxcr):
                maxcr = crlen

        # calculate jumps number and starting point
        jump_size = numpy.dot(norm_vec, v3)
        origin_shift_full = numpy.dot(self.origin, norm_vec)
        origin_shift = jump_size*math.modf(origin_shift_full/jump_size)[0]

        # for safety one jump before and after
        start_len = (math.ceil((mincr + origin_shift)/jump_size)-1)*jump_size
        jumps = math.floor((maxcr - start_len)/jump_size)+2
        start_vec = numpy.multiply(norm_vec, start_len)

        return (jumps, norm_vec, start_vec)


class Molecule:
    def __init__(self, types, coords):
        check_if_numpy2D(coords, "Molecule.__init__")
        self.coords = coords
        self.coords_n = len(coords)
        self.types = types

    def rotate(self, axis, theta, center=[0, 0, 0]):
        rotation_matrix = math_tools.rotation_matrix(axis, theta)
        rotation_center = numpy.array([center[0], center[1], center[2]],
                                      dtype=float)

        log.info("Rotating molecule around [%f, %f, %f] by %f radians"
                 % (axis[0], axis[1], axis[2], theta))

        coords = self.coords
        for i in xrange(len(coords)):
            atom_position = numpy.copy(coords[i])
            atom_position -= rotation_center
            atom_position = numpy.dot(rotation_matrix, atom_position)
            atom_position += rotation_center
            coords[i] = atom_position

    def center_com(self):
        # com - center of mass
        log.warning("Not implemented yet...")

    def center_coc(self):
        """
        coc - center of coordinates
        """
        coords = self.coords
        average = numpy.average(coords, axis=0)

        for i in xrange(self.coords_n):
            self.coords[i] = coords[i] - average

    def center_atom(self, atom_id=0):
        # center on atom, default first
        log.warning("Not implemented yet...")

    def get_molecule(self):
        return (self.types, self.coords)

class Crystal:
    def __init__(self, molecule, lattice, active=True):
        self.lattice = lattice
        self.molecule = molecule

        try:
            self.nodes = lattice.get_nodes()
            self.nodes_n = len(self.nodes)
        except:
            log.error("Something went wrong while getting nodes form lattice")

        if active:
            self.active_list = numpy.ones(self.nodes_n, dtype=bool)
        else:
            self.active_list = numpy.zeros(self.nodes_n, dtype=bool)

    def union(self, function):
        self.validate_function(function)
        active_list = self.active_list
        nd = self.nodes
        for i in xrange(len(active_list)):
            active_list[i] = active_list[i] | function(nd[i][0], nd[i][1],
                                                       nd[i][2])

    def intersect(self, function):
        self.validate_function(function)
        active_list = self.active_list
        nd = self.nodes
        for i in xrange(len(active_list)):
            active_list[i] = active_list[i] & function(nd[i][0], nd[i][1],
                                                       nd[i][2])

    def subtract(self, function):
        self.validate_function(function)
        active_list = self.active_list
        nd = self.nodes
        for i in xrange(len(active_list)):
            active_list[i] = active_list[i] & (not function(nd[i][0],
                                                            nd[i][1], nd[i][2]))

    def validate_function(self, function):
        try:
            function(0, 0, 0)
        except:
            log.error("Passed function is ill defined, should be:"
                      " bool function(x, y, z)")

    def get_crystal(self):
        nodes = self.nodes
        nodes_n = self.nodes_n
        active_list = self.active_list

        nodes_final_n = sum(active_list)

        mol_types, mol_coords = self.molecule.get_molecule()
        mol_n = len(mol_coords)

        crystal_n = nodes_final_n*mol_n

        crystal_coords = numpy.empty((crystal_n, 3), dtype=float)
        crystal_types = [None]*crystal_n

        index_current = 0
        for i in xrange(nodes_n):
            if not active_list[i]:
                continue

            for j in xrange(mol_n):
                crystal_coords[index_current] = nodes[i] + mol_coords[j]
                crystal_types[index_current] = mol_types[j]
                index_current += 1

        return (crystal_types, crystal_coords)

class System:
    def __init__(self):
        pass

class shape:
    @staticmethod
    def inside_box(xmin, xmax, ymin, ymax, zmin, zmax):
        return lambda x, y, z: x >= xmin and x <= xmax and y >= ymin and \
            y <= xmax and z >= zmin and z <= zmax

    @staticmethod
    def outside_box(xmin, xmax, ymin, ymax, zmin, zmax):
        return lambda x, y, z: not (x >= xmin and x <= xmax and y >= ymin \
                                    and y <= xmax and z >= zmin and z <= zmax)

    @staticmethod
    def inside_sphere(x0, y0, z0, r):
        return lambda x, y, z: (x-x0)**2 + (y-y0)**2 + (z-z0)**2 <= r**2

    @staticmethod
    def outside_sphere(x0, y0, z0, r):
        return lambda x, y, z: (x-x0)**2 + (y-y0)**2 + (z-z0)**2 > r**2


