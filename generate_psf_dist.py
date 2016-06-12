#!/usr/bin/env python3
from astropy.io import ascii
import matplotlib.pyplot as plt
import os

# Configuration
med_path = "./data/psf_median/"
img_path = "./plots/psf_median_dist/"
grid_x = 3
grid_y = 3
headers = ["x", "y", "flux", "bg", "bg_rms", 
           "mx", "my", "mxy", "a4", "a5", "a6"]

for h in headers:
    fig = plt.figure()
    for x in range(grid_x):
        for y in range(grid_y):
            median = ascii.read(med_path + "grid_" + str(x) + "_" + str(y) + "_" + h + ".psf")
            ax = plt.subplot2grid((grid_x+1, grid_y+1), (y, x))
            ax.hist(median[h], 50, normed=False, histtype='stepfilled', range=(min(median[h]),max(median[h])))
            ax.set_title("Grid: (" + str(x) + ":" + str(y) + 
                         ") Sample size: " + str(len(median[h])), fontsize=16)
            ax.set_xlabel(h, fontsize=16)
            ax.minorticks_on()
            ax.grid(True)
            plt.tick_params(which='both', width=2)
            plt.tick_params(which='major', length=5)
            plt.ticklabel_format(style='sci', scilimits=(-3,4))
    fig.set_dpi(300)
    fig.set_size_inches(20, 20)
    plt.tight_layout()
    os.makedirs(img_path, exist_ok=True)
    fig.savefig(img_path + h + ".pdf", bbox_inches='tight')