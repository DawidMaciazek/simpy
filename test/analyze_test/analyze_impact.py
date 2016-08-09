from simpy import analyze
import scipy.ndimage as ndimage
import numpy as np
import matplotlib.pyplot as plt


def test_area(frames):
    coord = frames[0]['coord']
    #plt.plot(coord.transpose()[0], coord.transpose()[2], "o")
    #plt.show()
    ct = coord.transpose()
    x = ct[0]
    z = ct[2]

    heatmap, xedges, yedges = np.histogram2d(x, z, bins=50)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    plt.clf()
    plt.imshow(heatmap, extent=extent)
    plt.show()


#test_area(frames)

def draw_image(frames):
    coord = frames[0]['coord']
    # lets bin this
    img = np.zeros((5,5), dtype=float)
    img[1][1] = 0.1
    img[2][1] = 0.11
    img[1][4] = 0.2
    img[0][1] = 0.27
    img[2][0] = 0.3


    plt.imshow(img)
    plt.show()


def bining_test(frames, frame_no):
    coord = frames[frame_no]['coord']
    coord = coord[coord[:,1] < 4]
    coord = coord[coord[:,1] > -4]

    tcoord = coord.transpose()
    x = tcoord[0]
    z = tcoord[2]
    #plt.plot(x, z, "o")

    binsize = 0.15
    offset = 0.1
    xmin = np.amin(x) - offset
    xmax = np.amax(x) + offset

    zmin = np.amin(z) - offset
    zmax = np.amax(z) + offset

    xbin_range = np.arange(xmin, xmax, binsize)
    xbin_id = np.digitize(x, xbin_range)
    xbins = len(xbin_range)

    zbin_range = np.arange(zmin, zmax, binsize)
    zbin_id = np.digitize(z, zbin_range)
    zbins = len(zbin_range)

    atoms_img = np.zeros((xbins, zbins), dtype=float)

    print xbin_id[:10]
    print zbin_id[:10]

    for i in xrange(len(x)):
        atoms_img[xbin_id[i]-1][zbin_id[i]-1] += 1.0

    atoms_img_conv = ndimage.gaussian_filter(atoms_img, sigma=12 , order=0)
    plt.imshow(atoms_img_conv)
    plt.show()


iname = "Ag_50.lammpstrj"
ifile = analyze.Traj(iname)
frames = ifile.read(20)

for i in range(20):
    bining_test(frames, i)
#draw_image(frames)

