#!/usr/bin/env python3
import focus_series as fs
from astropy.io import ascii
import numpy as np
import os
import math
import time
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import time
import argparse

# Configuration
tsi_path = "./data/tsi_tab/"
img_path = "./plots/tsi_dist/"

parser = argparse.ArgumentParser(description='Generate tsi histogram plots.')
parser.add_argument(
    '--tsi-path', 
    help='Path to the tsi meta-data files', 
    default=tsi_path
)
parser.add_argument(
    '--img-path', 
    help='Path to directory where the plots will be saved', 
    default=img_path
)
args = parser.parse_args()

tsi_path = args.tsi_path
img_path = args.img_path

start = time.time()
tsi = fs.read_all_tsi(tsi_path)

i = 0
_max = None
for k in tsi:
    if not fs.check_if_const(np.histogram(tsi[k], 50)[0]):
        fig = plt.figure()
        ax = fig.gca()
        ax.hist(tsi[k], 50, normed=False, histtype='stepfilled')
        ax.set_title("Sample size: " + str(len(tsi[k])), fontsize=12)
        ax.set_xlabel(k, fontsize=10)
        ax.minorticks_on()
        ax.grid(True)
        plt.tick_params(which='both', width=2)
        plt.tick_params(which='major', length=5)
        plt.ticklabel_format(style='sci', scilimits=(-3,4))
        fig.set_dpi(300)
        fig.set_size_inches(6, 6)
        plt.tight_layout()
        fig.savefig(img_path + k + ".pdf", dpi='figure', bbox_inches='tight')
        i += 1
    if (_max is not None and i >= _max):
        break
print(i, "Plots")
print(time.time()-start)