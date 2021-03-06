#!/usr/bin/env python3
import focus_series as fs
import argparse
import os

parser = argparse.ArgumentParser(
    description='Generate tsi heatmap of mutual information.'
)
parser.add_argument(
    '-f', '--file',
    help='File containing the mutual information',
    required=True
)
parser.add_argument(
    '-n', '--fields-per-plot',
    help='The number of parameters in each plot',
    default=20,
    type=int
)
parser.add_argument(
    '--plot-path',
    help='Path to the directory where the heatmap plots will be saved',
    default='./plots/'
)
parser.add_argument(
    '-s,', '--size',
    help='Size of the plot in inches. Default is automatic',
    default=None,
    type=float
)
parser.add_argument(
    '--type-1', 
    choices=['psf', 'tsi'], 
    help='Parameter type of the first key in the mutual info file',
    required=True
)
parser.add_argument(
    '--type-2',
    choices=['psf', 'tsi'],
    help='Parameter type of the second key in the mutual info file',
    required=True
)
args = parser.parse_args()
os.makedirs(args.plot_path, exist_ok=True)

if args.size != None:
    size = (args.size, args.size)
else:
    size = None

heatmap_kwargs = {}
if args.type_1 == "tsi":
    heatmap_kwargs['xlabel_callback'] = fs.get_sparse_tsi
if args.type_2 == "tsi":
    heatmap_kwargs['ylabel_callback'] = fs.get_sparse_tsi

fs.mutual_info_heatmap(args.file, args.plot_path, args.fields_per_plot, size,
                    **heatmap_kwargs)