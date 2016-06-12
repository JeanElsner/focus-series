#!/usr/bin/env python3
import focus_series as fs
import datetime as dt
import pytz
import time
import numpy as np
import json
import os

start = time.time()
tsi_path = './data/tsi_tab/'
meteo_path_src = './data/meteo_logs/'
meteo_path_dst = './data/meteo_preproc/'
time_delta = 150
tz_ger = pytz.timezone('Europe/Berlin')
os.makedirs(meteo_path_dst, exist_ok=True)
limit = 1000
i = 0
ff = fs.get_tsi_filenames(tsi_path)
print('Processing', len(ff), 'files')
for f in ff:
    mf = meteo_path_dst + f[:-4]
    if os.path.isfile(mf):
        print(f, 'already processed, skipping..')
        continue
    tsi = fs.read_tsi(tsi_path + f)
    tsi_dt = tz_ger.localize(dt.datetime.fromtimestamp(tsi['POSITION.LOCAL.UTC']))
    i += 1
    meteo = fs.search_meteo(tsi_dt, meteo_path_src, time_delta=time_delta)
    json.dump(meteo, open(mf, 'w'), indent=0)
    if i >= limit:
        break
print(time.time()-start)
