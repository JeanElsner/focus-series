#!/usr/bin/env python3
import focus_series as fs
from astropy.io import ascii
import numpy as np
import os
import math
import time
from scipy.stats import entropy

# Configuration
tsi_path = "./data/tsi_tab/"
med_path = "./data/psf_median/"
dat_path = "./data/"
grid_x = 3
grid_y = 3
headers = ["x", "y", "flux", "bg", "bg_rms", 
           "mx", "my", "mxy", "a4", "a5", "a6"]

start = time.time()
psf = fs.read_all_psf_grids(med_path, headers, grid_x, grid_y)
tsi = fs.read_all_tsi(tsi_path)
_max = None
i = 0
entr = []
#TODO check for focus[2].realpos -bestfocus

for k in psf:
    dist = np.histogram(k[0], 50, range=(min(k[0]), max(k[0])))
    if not fs.check_if_const(dist[0]):
        for k2 in tsi:
            dist2 = np.histogram(tsi[k2], 50)
            if not fs.check_if_const(dist2[0]):
                mi = fs.mutual_info(dist[0], dist2[0])
                i += 1
                entr.append([mi, k[1] + "(" + str(k[2]) + ":" + str(k[3]) + ")", k2 ])

    if (_max is not None and i >= _max):
        print("max reached")
        break

entr.sort(key=lambda x: -x[0])
ascii.write(np.array(entr), dat_path + "psf_tsi_mutual_info", names=['entropy', 'key1', 'key2'])

print(i, "iterations")
print(time.time()-start)