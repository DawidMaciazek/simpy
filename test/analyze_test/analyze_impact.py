from simpy import analyze
import numpy as np
import matplotlib.pyplot as plt

iname = "Ag_50.lammpstrj"
ifile = analyze.Traj(iname)
frames = ifile.read(4)


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

draw_image(frames)
