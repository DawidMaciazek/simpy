from simpy import analyze
from simpy import write
import numpy

from scipy.fftpack import ifftn
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import ticker

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
    write.xyz(output).write(['Si']*len(surf_array), surf_array)

def write_surf_matrix(surf_matrix, output):
    xlen = len(surf_matrix)
    ylen = len(surf_matrix[0])

    surf_array = numpy.empty((xlen*ylen, 3))
    for x_index in xrange(xlen):
        for y_index in xrange(ylen):
            surf_index = x_index*xlen+y_index
            surf_array[surf_index][0] = x_index
            surf_array[surf_index][1] = y_index
            surf_array[surf_index][2] = surf_matrix[x_index][y_index]

    write.xyz(output).write(['H']*len(surf_array), surf_array)


# box period
p_min = -59.641
p_max = 59.841
p = p_max - p_min
act_box=[[-59.641, 59.841], [-59.641, 59.841]]

# angle = 85
surface_file = "../resources/85_90bond_last.xyz"
small_coords, types = get_atoms(surface_file)

repeat_units = [[-p, -p], [-p, 0.0], [-p, p],
                [0.0, -p], [0.0, 0.0], [0.0, p],
                [p, -p], [p, 0.0], [p, p]]

# repeat coords on 9x9 grig
coords =  numpy.empty((0, 3), dtype=float)
for repeat in repeat_units:
    coords = numpy.append(coords, small_coords+[0.0, repeat[0], repeat[1]], axis=0)

# Initialize class
# 1.4 - surface atom radius. 1.6 - probe atom radius,
rms = analyze.RMS(coords, 1.4, 1.6, sampling_box=act_box, offset=0.0)


# Set sampling array with spacing 1A
rms.calc_sampling_array(1.0)
surf_matrix = rms.get_surf_array_oa2d(3, format='matrix')

# Compute Fourier transform, then extract real and absolut part
fft_matrix = ifftn(surf_matrix)
fft_matrix_abs = numpy.absolute(fft_matrix)
fft_matrix_real = numpy.real(fft_matrix)

# Save
write_surf_matrix(fft_matrix_real, 'fft_85a_real.xyz')
write_surf_matrix(fft_matrix_abs, 'fft_85a_abs.xyz')
write_surf_matrix(surf_matrix, 'surf_85a.xyz')

# FT value range on images
fft_range = 0.05


# Visualize results
plt.figure(facecolor='white')
surf_zero = surf_matrix-surf_matrix.min()
sp = plt.subplot(2,3,1)
sp.set_ylabel(r'$\alpha = 85^o$', fontsize=26)
sp.set_title('Surface topology', fontsize=20)
im = plt.imshow(surf_zero, cmap=cm.afmhot)
cb = plt.colorbar(im, fraction=0.046, pad=0.04)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()


sp = plt.subplot(2,3,2)
sp.set_title('Fourier transform (Real part)', fontsize=20)
im = plt.imshow(fft_matrix_real, cmap=cm.afmhot, clim=(-fft_range, fft_range))
cb = plt.colorbar(im, fraction=0.046, pad=0.04)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

sp = plt.subplot(2,3,3)
sp.set_title('Fourier transform (Abs)', fontsize=20)
im = plt.imshow(fft_matrix_abs, cmap=cm.afmhot, clim=(-fft_range, fft_range))
cb = plt.colorbar(im, fraction=0.046, pad=0.04)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()


# ------------------------------------------------------------------
# angle = 89
surface_file = "../resources/89_90bond_last.xyz"
small_coords, types = get_atoms(surface_file)

repeat_units = [[-p, -p], [-p, 0.0], [-p, p],
                [0.0, -p], [0.0, 0.0], [0.0, p],
                [p, -p], [p, 0.0], [p, p]]

# repeat coords on 9x9 grig
coords =  numpy.empty((0, 3), dtype=float)
for repeat in repeat_units:
    coords = numpy.append(coords, small_coords+[0.0, repeat[0], repeat[1]], axis=0)

# Initialize class
# 1.4 - surface atom radius. 1.6 - probe atom radius,
rms = analyze.RMS(coords, 1.4, 1.6, sampling_box=act_box, offset=0.0)


rms.calc_sampling_array(1.0)
surf_matrix = rms.get_surf_array_oa2d(3, format='matrix')


fft_matrix = ifftn(surf_matrix)
fft_matrix_abs = numpy.absolute(fft_matrix)
fft_matrix_real = numpy.real(fft_matrix)

write_surf_matrix(fft_matrix_real, 'fft_89a_real.xyz')
write_surf_matrix(fft_matrix_abs, 'fft_89a_abs.xyz')
write_surf_matrix(surf_matrix , 'surf_89a.xyz')



surf_zero = surf_matrix-surf_matrix.min()
sp = plt.subplot(2,3,4)
sp.set_ylabel(r'$\alpha = 89^o$', fontsize=26)
sp.set_title('Surface topology', fontsize=20)
im = plt.imshow(surf_zero, cmap=cm.afmhot)
cb = plt.colorbar(im, fraction=0.046, pad=0.04)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()


sp = plt.subplot(2,3,5)
sp.set_title('Fourier transform (Real part)', fontsize=20)
im = plt.imshow(fft_matrix_real, cmap=cm.afmhot, clim=(-fft_range, fft_range))
cb = plt.colorbar(im, fraction=0.046, pad=0.04)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()

sp = plt.subplot(2,3,6)
sp.set_title('Fourier transform (Abs)', fontsize=20)
im = plt.imshow(fft_matrix_abs, cmap=cm.afmhot, clim=(-fft_range, fft_range))
cb = plt.colorbar(im, fraction=0.046, pad=0.04)
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()


plt.tight_layout()
plt.show()
