from simpy import analyze
from scipy.signal import find_peaks_cwt
from scipy.optimize import curve_fit
import numpy as np

def bin_coord(coord, lo, hi, bin_size):
    bin_number = (hi-lo)/float(bin_size) + 1

    bin_centers = np.linspace(lo, hi, bin_number)
    bin_edges = bin_centers[:-1] + 0.5*bin_size

    bin_indexes = np.digitize(coord, bin_edges)

    bin_cnt = np.zeros(bin_centers.shape, dtype=int)
    for i in bin_indexes:
        bin_cnt[i] += 1
    #dict(zip(bin_centers, bin_cnt))
    return bin_centers, bin_cnt

def sigmoid(x, x0, k, a, b):
    y = a / (1.0 + np.exp(-k*(x-x0))) + b
    return y

traj = analyze.Traj("big_05501.lammpstrj")
#traj = analyze.Traj("85deg.lammpstrj")
fr = traj.read()

cr = fr["coord"]
bin_centers, bin_cnt = bin_coord(cr[:,2], -160, 20, 0.1)
print(bin_centers.shape, bin_cnt.shape)
bin_cnt = (bin_cnt + np.roll(bin_cnt, 1) + np.roll(bin_cnt, -1))/3.0
save_data = np.transpose(np.array([bin_centers, bin_cnt]))

np.savetxt("out.txt", save_data)
peaks = find_peaks_cwt(bin_cnt, np.arange(5,10))
np.savetxt("peaks.txt", save_data[peaks])
#binned = dict(zip(bin_centers, bin_cnt))

xdata = bin_centers[peaks]
ydata = bin_cnt[peaks]

init_optimal = [-30,-0.1,300,100]
coef, pcov = curve_fit(sigmoid, xdata, ydata, p0=init_optimal)

fit_curve = sigmoid(xdata, coef[0], coef[1], coef[2], coef[3])
print(coef)

np.savetxt("fit.txt",  np.transpose(np.array([xdata, fit_curve])))
