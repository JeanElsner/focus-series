#!/usr/bin/env python3
import focus_series as fs
import numpy as np
import json
import matplotlib.pyplot as plt
import time
import math
import argparse
from astropy.io import ascii
from scipy import interpolate

psf_path   = './data/psf_data/'
meteo_path = './data/meteo_preproc/'
plot_path  = './plots/wind/'

velocity    = 'WIND.VELOCITY.1h.AVG'
direction   = 'WIND.DIRECTION.1h.AVG'
temperature = 'TEMPERATURE.5m.AVG'

parser = argparse.ArgumentParser(description='Analyse influence of wind on psf.')
parser.add_argument(
    '-p', '--parameter', 
    help='The psf parameter to analyse', 
    default='a4'
)
parser.add_argument(
    '--psf-path', 
    help='Path to the psf files', 
    default=psf_path
)
parser.add_argument(
    '--meteo_path', 
    help='Path to the preprocessed meteological data files', 
    default=meteo_path
)
parser.add_argument(
    '--plot_path',
    help='Path to directory where the plots will be saved',
    default=plot_path
)
parser.add_argument(
    '--size-multiplier',
    help='Multiplier for the size of the psf scatter plot',
    default=10,
    type=float
)
args = parser.parse_args()

start = time.time()
theta, r, a4, temp = ([],[],[],[])

for f in fs.get_psf_filenames(args.psf_path):
    psf = ascii.read(args.psf_path + f)
    try:
        meteo = json.load(open(args.meteo_path + f[:-4], 'r'))
    except FileNotFoundError:
        continue
    try:
        wind_velo = float(meteo[velocity].replace('m/s', '').strip())
        wind_dir = float(meteo[direction].replace('deg', '').strip())
        theta.append(np.radians(wind_dir))
        r.append(wind_velo)
        a4.append(np.nanmean(np.array(psf[args.parameter], dtype=float)))
        temp.append(float(meteo[temperature].replace('degC', '').strip()))
    except KeyError:
        print('No wind measurements for', f)
print(len(a4), 'datapoints')

_a4, knots    = fs.reduce_to_grid(theta, r, a4, 12, 12, lambda z,x,y: float('nan'))
_temp, knots2 = fs.reduce_to_grid(theta, r, temp, 12, 12, lambda x,y,z: float('nan'))

points = np.meshgrid(knots[0], knots[1], indexing='ij')
points2 = np.meshgrid(knots2[0], knots2[1], indexing='ij')

grid_x, grid_y = np.mgrid[0:2*3.14:100j, 0:20:100j]

""" Note: interpolation doesn't work at this point
    without resorting to extrapolation (cf. fallback in
    reduce_to_grid). A proposed solution would be to 
    interpolate a grid in polar coordinates before
    coordinate transformation.
"""
spline = interpolate.SmoothBivariateSpline(
    np.array(points[0]).flatten(), 
    np.array(points[1]).flatten(), 
    np.array(_a4).flatten(), 
    kx=3, ky=3
)
grid_z = spline.ev(grid_x, grid_y)

spline2 = interpolate.SmoothBivariateSpline(
    np.array(points2[0]).flatten(), 
    np.array(points2[1]).flatten(), 
    np.array(_temp).flatten(), 
    kx=3, ky=3
)
grid_z2 = spline2.ev(grid_x, grid_y)


fig = plt.figure()
fig.set_size_inches(14, 14)

# Contour plot psf
ax = plt.subplot(221, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title(args.parameter + ' contour', y=1.08)
ax.set_xlabel('south')
ax.text(1.1, 0.5, 'east', transform=ax.transAxes, verticalalignment='center')
plt.contourf(grid_x, grid_y, grid_z, alpha=.85)
cs = plt.contourf(points[0], points[1], abs(np.array(_a4)), alpha=.85)
cbar = plt.colorbar(cs, shrink=.8, pad=.15)
cbar.ax.set_ylabel(args.parameter)

# Scatter plot psf
ax = plt.subplot(223, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title(args.parameter + ' scatter, ' + str(len(a4)) + ' samples', y=1.08)
ax.set_xlabel('south')
ax.text(1.1, 0.5, 'east', transform=ax.transAxes, verticalalignment='center')
cs2 = plt.scatter(theta, r, c=abs(np.array(a4)), s=abs(np.array(a4)*args.size_multiplier), cmap=plt.cm.coolwarm, alpha=.75)
cbar2 = plt.colorbar(cs2, shrink=.8, pad=.15)
cbar2.ax.set_ylabel(args.parameter)

# Contour plot temperature
ax = plt.subplot(222, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title('Temperature contour', y=1.08)
ax.set_xlabel('south')
ax.text(1.1, 0.5, 'east', transform=ax.transAxes, verticalalignment='center')
plt.contourf(grid_x, grid_y, grid_z2, alpha=.85)
cs3 = plt.contourf(points2[0], points2[1], _temp, alpha=.85)
cbar3 = plt.colorbar(cs3, shrink=.8, pad=.15)
cbar3.ax.set_ylabel('Temperature')

# Scatter plot temperature
ax = plt.subplot(224, projection='polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_title('Temperature scatter, ' + str(len(temp)) + ' sampels', y=1.08)
ax.set_xlabel('south')
ax.text(1.1, 0.5, 'east', transform=ax.transAxes, verticalalignment='center')
cs4 = plt.scatter(theta, r, c=np.array(temp), s=np.array(temp)*5, cmap=plt.cm.coolwarm, alpha=.75)
cbar4 = plt.colorbar(cs4, shrink=.8, pad=.15)
cbar4.ax.set_ylabel('Temperature')

plt.tight_layout()
plt.savefig(args.plot_path + 'wind_' + args.parameter + '.pdf')
print(time.time()-start, 'seconds')
