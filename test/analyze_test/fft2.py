from simpy import analyze
from simpy import write
import numpy


def get_atoms(xyz_file):
    frames = analyze.Traj(xyz_file, 'xyz').read(1)

    coords = analyze.get(frames[0], 'coords')
    types = analyze.get(frames[0], 'element')

    return coords, types


def check_parameters(rms):
    print "active box size after auto calc:"
    print str(rms.active_box)

    print "\nr_sample, r_probe, r_vdw"
    print rms.r_sample, rms.r_probe, rms.r

    print "\noffset:"
    print str(rms.offset)


def write_surf_array(surf_array, output):
    write.xyz(output).write(['H']*len(surf_array), surf_array)

def write_surf_matrix(surf_matrix, output):
    xlen = len(surf_matrix)
    ylen = len(surf_matrix[0])

    surf_array = numpy.empty((xlen*ylen, 3))
    for i in xrange(xlen):
        for j in xrange(ylen):
            surf_index = i*xlen+j
            surf_array[surf_index][0] = i
            surf_array[surf_index][1] = j
            surf_array[surf_index][2] = surf_matrix[i][j]

    write.xyz(output).write(['He']*len(surf_array), surf_array)


surface_file = "../resources/surf_test.xyz"
coords, types = get_atoms(surface_file)

rms = analyze.RMS(coords, 1.4, 1.6)
rms.calc_sampling_array(1.0)

surf_matrix = rms.get_surf_array_oa2d(3, format='matrix')


# plot results

from scipy.fftpack import ifftn
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import gridspec

fft_matrix = numpy.real(ifftn(surf_matrix))

fig = plt.figure(figsize=(8, 4))
grid = gridspec.GridSpec(1, 2, width_ratios=[1, 1])
ax0 = plt.subplot(grid[0])
ax0.set_title("Surface topology")
ax0.imshow(surf_matrix, cmap=cm.afmhot)

ax1 = plt.subplot(grid[1])
ax1.set_title("Fourier transform (Real part)")
ax1.imshow(fft_matrix, cmap=cm.afmhot)


plt.tight_layout()
plt.show()
