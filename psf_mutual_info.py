#!/usr/bin/env python3
import focus_series as fs
from astropy.io import ascii
import numpy as np
import time
from scipy.stats import entropy

# Configuration
med_path = "./data/psf_median/"
dat_path = "./data/"
grid_x = 3
grid_y = 3
headers = ["x", "y", "flux", "bg", "bg_rms", 
           "mx", "my", "mxy", "a4", "a5", "a6"]

start = time.time()
psf = fs.read_all_psf_grids(med_path, headers, grid_x, grid_y)
entr = []
_max = None
i = 0

for k in psf:
    dist = np.histogram(k[0], 50, range=(min(k[0]), max(k[0])))
    if not fs.check_if_const(dist[0]):
        for k2 in psf:
            dist2 = np.histogram(k2[0], 50, range=(min(k2[0]), max(k2[0])))
            if not fs.check_if_const(dist2[0]):
                mi = fs.mutual_info(dist[0], dist2[0])
                if k[1] != k2[1]:
                    i += 1
                    entr.append([mi, k[1] + "(" + str(k[2]) + ":" + str(k[3]) + ")", 
                                 k2[1] + "(" + str(k2[2]) + ":" + str(k2[3]) + ")"])

    if (_max is not None and i >= _max):
        print("maximum reached")
        break

entr.sort(key=lambda x: -x[0])
ascii.write(np.array(entr), dat_path + "psf_mutual_info", names=['mutual_info', 'psf1', 'psf2'])

print(i, "iterations")
print(time.time()-start)