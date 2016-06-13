#!/usr/bin/env python3
import focus_series as fs
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import os
from astropy.io import ascii
sns.set(context='paper', font='monospace')
plt.ioff()

parser = argparse.ArgumentParser(
    description='Generate tsi heatmap of mutual information.'
)
parser.add_argument(
    '-f', '--file',
    help='File containing the mutual information',
    default='./data/tsi_mutual_info'
)
parser.add_argument(
    '-n', '--fields-per-plot',
    help='The number of tsi fields in each plot',
    default=20,
    type=int
)
parser.add_argument(
    '--plot-path',
    help='Path to the directory where the heatmap plots will be saved',
    default='./plots/tsi_heatmaps/'
)
parser.add_argument(
    '-s,', '--size',
    help='Size of the plot in inches. Default is automatic',
    default=None,
    type=float
)
args = parser.parse_args()
os.makedirs(args.plot_path, exist_ok=True)

if args.size != None:
    size = (args.size, args.size)
else:
    size = None

fs.mutual_info_heatmap(args.file, args.plot_path, args.fields_per_plot, size,
                    key1='tsi1', key2='tsi2',
                    xlabel_callback=fs.get_sparse_tsi,
                    ylabel_callback=fs.get_sparse_tsi)