import numpy
import math
import random


def rotation_matrix(r_vec, ra):
    r_vec = numpy.array(r_vec)
    r_vec = r_vec/numpy.linalg.norm(r_vec)
    rmatrix = [[math.cos(ra) + r_vec[0]*r_vec[0]*(1-math.cos(ra)),
                r_vec[0]*r_vec[1]*(1-math.cos(ra)) - r_vec[2]*math.sin(ra),
                r_vec[0]*r_vec[2]*(1-math.cos(ra)) + r_vec[1]*math.sin(ra)],
               [r_vec[1]*r_vec[0]*(1-math.cos(ra)) + r_vec[2]*math.sin(ra),
                math.cos(ra) + r_vec[1]*r_vec[1]*(1-math.cos(ra)),
                r_vec[1]*r_vec[2]*(1-math.cos(ra)) - r_vec[0]*math.sin(ra)],
               [r_vec[2]*r_vec[0]*(1-math.cos(ra)) - r_vec[1]*math.sin(ra),
                r_vec[2]*r_vec[1]*(1-math.cos(ra)) + r_vec[0]*math.sin(ra),
                math.cos(ra) + r_vec[2]*r_vec[2]*(1-math.cos(ra))]]

    return numpy.array(rmatrix)


def random_unit_sphere(numpy_arr=True):
        angle = random.random()*math.pi*2
        u = random.uniform(-1, 1)
        u2 = u*u

        x = math.sqrt(1-u2) * math.cos(angle)
        y = math.sqrt(1-u2) * math.sin(angle)
        z = u

        if numpy_arr:
            return numpy.array([x, y, z], dtype=float)
        else:
            return [x, y, z]
