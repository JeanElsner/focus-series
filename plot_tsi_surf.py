#!/usr/bin/env python3
import focus_series as fs
import matplotlib.pyplot as plt
import numpy as np
import time
import os
import argparse
from matplotlib import cm
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D

# Configuration
tsi_path = "./data/tsi_tab/"
img_path = "./plots/tsi_surface/"

VAL = 'POSITION.INSTRUMENTAL.FOCUS[2].REALPOS'

grid_temp = 4
grid_elev = 4

kind_dict = {'std' : "Standard Deviation", 'mean' : "Mean", 'med' : "Median", 'var' : "Variance"}
var_dict  = {'POSITION.INSTRUMENTAL.FOCUS[2].REALPOS' : "FOCUS[2].REALPOS"}

parser = argparse.ArgumentParser(description='Calculate tsi surface interpolant.')
parser.add_argument(
    '-k', '--kind', 
    help='The kind of average to calculate (default: mean)',
    choices=['mean', 'std', 'var', 'med'],
    default='mean'
)
args = parser.parse_args()
kind = args.kind

start = time.time()

temp_elev_lst = fs.read_all_tsi_val_temp_elev(VAL, tsi_path)
gridded = fs.populate_temp_elev_dict(temp_elev_lst, grid_temp, grid_elev, single_value=True)

fig = plt.figure()
fig.set_size_inches(14, 14)

grid_x, grid_y = np.mgrid[min(gridded['temp']):max(gridded['temp']):200j, 
                          min(gridded['elev']):max(gridded['elev']):200j]
spline = interpolate.SmoothBivariateSpline(gridded['temp'], gridded['elev'], gridded[kind])

print('knots\n', spline.get_knots())
print('coefficients\n', spline.get_coeffs())

grid_z = spline.ev(grid_x, grid_y)
i = len(temp_elev_lst['val'])

"""
grid_z2 = []
for u in np.linspace(min(gridded['temp']), max(gridded['temp']), 200):
    for v in np.linspace(min(gridded['elev']), max(gridded['elev']), 200):
        grid_z2.append(eval_surf_spline(u, v, knots[0], knots[1], coeffs))
        
grid_z2 = np.reshape(grid_z2, (200, 200))
grid_x2, grid_y2 = np.mgrid[min(gridded['temp']):max(gridded['temp']):200j, 
                          min(gridded['elev']):max(gridded['elev']):200j]
#print(grid_z2)
"""
        
# First plot
ax = fig.add_subplot(221, projection='3d')
ax.set_ylabel("ZD")
ax.set_zlabel(var_dict[VAL])
ax.set_title(var_dict[VAL] + " " + kind_dict[kind] + ", " + str(i) + " samples")
ax.azim = 0
ax.elev = 0

ax.plot_surface(grid_x, grid_y, grid_z, zorder=10, cmap=cm.coolwarm, antialiased=True, shade=True, alpha=.25)
ax.scatter(gridded['temp'], gridded['elev'], gridded[kind], c='r', zorder=1)
ax.scatter(temp_elev_lst['temp'], temp_elev_lst['elev'], temp_elev_lst['val'], c='b', zorder=1)

# Second plot
ax = fig.add_subplot(222, projection='3d')
ax.set_xlabel("Temperature")
ax.set_zlabel(var_dict[VAL])
ax.set_title(var_dict[VAL] + " " + kind_dict[kind] + ", " + str(i) + " samples")
ax.azim = -90
ax.elev = 0

ax.plot_surface(grid_x, grid_y, grid_z, zorder=10, cmap=cm.coolwarm, antialiased=True, shade=True, alpha=.25)
ax.scatter(gridded['temp'], gridded['elev'], gridded[kind], c='r', zorder=1)
ax.scatter(temp_elev_lst['temp'], temp_elev_lst['elev'], temp_elev_lst['val'], c='b', zorder=1)

# Third Plot
ax = fig.add_subplot(223, projection='3d')
ax.set_xlabel("Temperature")
ax.set_ylabel("ZD")
ax.set_zlabel(var_dict[VAL])
ax.set_title(var_dict[VAL] + " " + kind_dict[kind] + ", " + str(i) + " samples")
ax.azim = -60
ax.elev = 20

#ax.scatter(grid_x2, grid_y2, grid_z2, c='r', zorder=1)
#ax.plot_surface(grid_x2, grid_y2, grid_z2, zorder=10, cmap=cm.coolwarm, antialiased=True, shade=True, alpha=.75)
ax.plot_surface(grid_x, grid_y, grid_z, zorder=10, cmap=cm.coolwarm, antialiased=True, shade=True, alpha=.75)
#ax.scatter(gridded['temp'], gridded['elev'], gridded[kind], c='r', zorder=1)
#ax.scatter(temp_elev_lst['temp'], temp_elev_lst['elev'], temp_elev_lst['val'], c='b', zorder=1)
#cset = ax.contourf(grid_x, grid_y, grid_z, cmap=cm.coolwarm, zdir='z', offset=min(temp_elev_lst['val']), zorder=2)

# Fourth Plot
ax = fig.add_subplot(224, projection='3d')
ax.set_xlabel("Temperature")
ax.set_ylabel("ZD")
ax.set_zlabel(var_dict[VAL])
ax.set_title(var_dict[VAL] + " " + kind_dict[kind] + ", " + str(i) + " samples")
ax.azim = -60
ax.elev = 20

from matplotlib import cm
ax.plot_surface(grid_x, grid_y, grid_z, zorder=10, cmap=cm.coolwarm, antialiased=True, shade=True, alpha=.75)
ax.scatter(gridded['temp'], gridded['elev'], gridded[kind], c='r', zorder=1)

for i in range(len(gridded['err'])):
    ax.plot([gridded['temp'][i], gridded['temp'][i]], [gridded['elev'][i], gridded['elev'][i]], 
            [gridded[kind][i]-gridded['err'][i], gridded[kind][i]+gridded['err'][i]], 
            marker='_', linewidth=2.0, color='red')

plt.tight_layout()
os.makedirs(img_path, exist_ok=True)
fig.savefig(img_path + VAL + "_" + kind + ".pdf", dpi='figure', bbox_inches='tight')

print(time.time()-start)