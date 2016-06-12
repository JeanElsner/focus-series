#!/usr/bin/env python3
from astropy.io import ascii
import focus_series as fs
import numpy as np
import os
import math
import time

# Configuration
psf_path = "./data/psf_data/"
med_path = "./data/psf_median/"
max_x = 9000
max_y = 9000
grid_x = 3
grid_y = 3
headers = ["x", "y", "flux", "bg", "bg_rms", 
           "mx", "my", "mxy", "a4", "a5", "a6"]

start = time.time()
os.makedirs(med_path, exist_ok=True)
median = fs.populate_psf_grid(psf_path, grid_x, grid_y, max_x, max_y, headers)

for x in range(grid_x):
    for y in range(grid_y):
        for c in median[x][y]:
            ascii.write([median[x][y][c]], med_path + "grid_" + str(x) 
                        + "_" + str(y) + "_" + c + ".psf", names=[c])

print(time.time()-start)