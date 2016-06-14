#!/usr/bin/env python3
import focus_series as fs
import numpy as np
import os
import time
from astropy.io import ascii

# Configuration
tsi_path = "./data/tsi_tab/"
dat_path = "./data/"

start = time.time()

tsi = fs.read_all_tsi(tsi_path)
entr = []
_max = None
i = 0

for k in tsi:
    dist = np.histogram(tsi[k], 50)
    if not fs.check_if_const(dist[0]):
        for k2 in tsi:
            dist2 = np.histogram(tsi[k2], 50)
            if not fs.check_if_const(dist2[0]):
                mi = fs.mutual_info(dist[0], dist2[0])
                if True:
                #if k != k2:
                    i += 1
                    entr.append([mi, k, k2])
    if (_max is not None and i >= _max):
        print("max reached")
        break

entr.sort(key=lambda x: -x[0])
os.makedirs(dat_path, exist_ok=True)
ascii.write(np.array(entr), dat_path + "tsi_mutual_info", names=['mutual_info', 'key1', 'key2'])

print(i, "iterations")
print(time.time()-start)