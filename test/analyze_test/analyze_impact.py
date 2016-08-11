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

    element = frames[frame_no]['element']
    selel = np.array([True if el == 'Ag' else False for el in element])
    strip = 4
    #xmax = 45
    xmax = 100
    zmax = 5
    coord = frames[frame_no]['coord']
    sely = np.logical_and(coord[:,1] < strip, coord[:,1] > -strip)
    selx = np.logical_and(coord[:,0] < xmax, coord[:,0] > -xmax)
    selz = np.logical_and(coord[:,2] < zmax, coord[:,2] > -xmax)
    sel = np.logical_and(np.logical_and(selx, sely), selz)
    sel = np.logical_and(sel, selel)


    ek = frames[frame_no]['c_sput_ke']

    element = np.array(element)
    he_sel = np.logical_and(sel, ek > 18.0)
    high_ek = element[he_sel]

    ek = ek[sel]
    print high_ek
    print coord[he_sel]

    ep = frames[frame_no]['c_sput_pe']
    ep = ep[sel]

    coord = coord[sel]
    ocoord = frames[0]['coord'][sel]

    tcoord = coord.transpose()
    x = tcoord[0]
    z = tcoord[2]
    #plt.plot(x, z, "o")

    binsize = 0.25
    offset = 0.1
    xmin = -xmax  - offset
    xmax = xmax + offset

    zmin = -xmax - offset
    zmax = zmax + offset

    xbin_range = np.arange(xmin, xmax, binsize)
    xbin_id = np.digitize(x, xbin_range)
    xbins = len(xbin_range)

    zbin_range = np.arange(zmin, zmax, binsize)
    zbin_id = np.digitize(z, zbin_range)
    zbins = len(zbin_range)

    atoms_img = np.zeros((xbins, zbins), dtype=float)

    energy_flag = False
    occupy_flag = True
    if energy_flag:
        e = ep
        for i in xrange(len(x)):
            e_loc = e[i]+0.11
            atoms_img[xbin_id[i]-1][zbin_id[i]-1] += e_loc
        print '[', np.amax(e), ',',  np.amin(e), ']'
        print np.mean(e), '(', np.std(e), ')'
    elif occupy_flag:
        for i in xrange(len(x)):
            atoms_img[xbin_id[i]-1][zbin_id[i]-1] = 1.0
        print np.amax(atoms_img)
    else:
        max_d = 0.0
        for i in xrange(len(x)):
            d = np.linalg.norm(coord[i]-ocoord[i])
            if d > max_d:
                max_d = d
            if d > 5:
                d = 5
            atoms_img[xbin_id[i]-1][zbin_id[i]-1] += d
        print max_d

    #plt.hist(e, 50, range=maxy)
    #plt.show()

    atoms_img_conv = ndimage.gaussian_filter(atoms_img, sigma=6 , order=0)
    plt.imshow(atoms_img_conv)
    plt.show()




iname = "Ag_50.lammpstrj"
iname = "/home/dawid/work/improvingSim/wave/simpleTest/Ag_100.lammpstrj"
ifile = analyze.Traj(iname)
frames = ifile.read(40)


print frames[0].keys()
for i in range(len(frames)):
    bining_test(frames, i)
#draw_image(frames)

