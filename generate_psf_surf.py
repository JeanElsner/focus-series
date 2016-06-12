#!/usr/bin/env python3
import focus_series as fs
import numpy as np
import os
import time
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from scipy import interpolate
import argparse

# Configuration
psf_path = "./data/psf_data/"
tsi_path = "./data/tsi_tab/"
img_path = "./plots/psf_surface/"
headers = ["flux", "bg", "bg_rms", 
           "mx", "my", "mxy", "a4", "a5", "a6"]

kind_dict  = {'std' : "Standard Deviation", 'mean' : "Mean", 'med' : "Median", 'var' : "Variance"}
headers    = ["a4", "a5", "a6"]
grid_temp  = 4
grid_elev  = 4

parser = argparse.ArgumentParser(description='Calculate interpolant psf surface.')
parser.add_argument(
    '-k', '--kind', 
    help='The kind of average to calculate (default: mean)',
    choices=['mean', 'std', 'var', 'med'],
    default='mean'
)
args = parser.parse_args()
kind = args.kind

start = time.time()
temp_elev_lst = fs.read_all_psf_temp_elev(headers, psf_path, tsi_path)
gridded = fs.generate_temp_elev_header_dict(headers)

for h in headers:
    gridded[h] = fs.populate_temp_elev_dict(temp_elev_lst[h], grid_temp, grid_elev, single_value=False)

for h in headers:
    i = len(temp_elev_lst[h][kind])
    fig = plt.figure()
    fig.set_size_inches(14, 14)
    
    grid_x, grid_y = np.mgrid[min(gridded[h]['temp']):max(gridded[h]['temp']):200j, 
                              min(gridded[h]['elev']):max(gridded[h]['elev']):200j]
    spline = interpolate.SmoothBivariateSpline(gridded[h]['temp'], gridded[h]['elev'], gridded[h][kind])
    grid_z = spline.ev(grid_x, grid_y)
    
    # First plot
    ax = fig.add_subplot(221, projection='3d')
    ax.set_ylabel("ZD")
    ax.set_zlabel(h)
    ax.set_title(h + " " + kind_dict[kind] + ", " + str(i) + " samples")
    ax.azim = 0
    ax.elev = 0
    
    ax.plot_surface(grid_x, grid_y, grid_z, zorder=10, cmap=cm.coolwarm, antialiased=True, shade=True, alpha=.25)
    ax.scatter(gridded[h]['temp'], gridded[h]['elev'], gridded[h][kind], c='r', zorder=1)
    ax.scatter(temp_elev_lst[h]['temp'], temp_elev_lst[h]['elev'], temp_elev_lst[h][kind], c='b', zorder=1)

    # Second plot
    ax = fig.add_subplot(222, projection='3d')
    ax.set_xlabel("Temperature")
    ax.set_zlabel(h)
    ax.set_title(h + " " + kind_dict[kind] + ", " + str(i) + " samples")
    ax.azim = -90
    ax.elev = 0
    
    ax.plot_surface(grid_x, grid_y, grid_z, zorder=10, cmap=cm.coolwarm, antialiased=True, shade=True, alpha=.25)
    ax.scatter(gridded[h]['temp'], gridded[h]['elev'], gridded[h][kind], c='r', zorder=1)
    ax.scatter(temp_elev_lst[h]['temp'], temp_elev_lst[h]['elev'], temp_elev_lst[h][kind], c='b', zorder=1)

    # Third Plot
    ax = fig.add_subplot(223, projection='3d')
    ax.set_xlabel("Temperature")
    ax.set_ylabel("ZD")
    ax.set_zlabel(h)
    ax.set_title(h + " " + kind_dict[kind] + ", " + str(i) + " samples")
    ax.azim = -60
    ax.elev = 20
    
    ax.plot_surface(grid_x, grid_y, grid_z, zorder=10, cmap=cm.coolwarm, antialiased=True, shade=True, alpha=.75)
    ax.scatter(gridded[h]['temp'], gridded[h]['elev'], gridded[h][kind], c='r', zorder=1)
    ax.scatter(temp_elev_lst[h]['temp'], temp_elev_lst[h]['elev'], temp_elev_lst[h][kind], c='b', zorder=1)
    
    try:
        from matplotlib import cm
        cset = ax.contourf(grid_x, grid_y, grid_z, cmap=cm.coolwarm, zdir='z', offset=min(temp_elev_lst[h][kind]), zorder=2)
    except:
        print('warning')
    
    # Fourth Plot
    ax = fig.add_subplot(224, projection='3d')
    ax.set_xlabel("Temperature")
    ax.set_ylabel("ZD")
    ax.set_zlabel(h)
    ax.set_title(h + " " + kind_dict[kind] + ", " + str(i) + " samples")
    ax.azim = -140
    ax.elev = 20
    
    ax.plot_surface(grid_x, grid_y, grid_z, zorder=10, cmap=cm.coolwarm, antialiased=True, shade=True, alpha=.75)
    ax.scatter(gridded[h]['temp'], gridded[h]['elev'], gridded[h][kind], c='r', zorder=1)
    
    for j in range(len(gridded[h]['err'])):
        ax.plot([gridded[h]['temp'][j], gridded[h]['temp'][j]], [gridded[h]['elev'][j], gridded[h]['elev'][j]], 
                [gridded[h][kind][j]-gridded[h]['err'][j], gridded[h][kind][j]+gridded[h]['err'][j]], 
                marker='_', linewidth=2.0, color='red')
    
    plt.tight_layout()
    os.makedirs(img_path, exist_ok=True)
    fig.savefig(img_path + h + "_" + kind + ".pdf", dpi='figure', bbox_inches='tight')

print(time.time()-start)