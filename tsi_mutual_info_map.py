#!/usr/bin/env python3
import focus_series as fs
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import pandas
import math
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

mi = ascii.read(args.file)
fields = sorted(set(list(mi['tsi1']) + list(mi['tsi2'])))
f_len  = len(fields)
f_num  = args.fields_per_plot

print(f_len, 'parameters total')
for i in range(math.ceil(f_len/f_num)):
    start = i*f_num
    end   = min((i+1)*f_num,f_len)
    sub_fields = fields[start:end]

    for j in range(math.ceil(f_len/f_num)):
        start2 = j*f_num
        end2   = min((j+1)*f_num,f_len)
        sub_fields2 = fields[start2:end2]
        mat = {k : {k : 0 for k in sub_fields2} for k in sub_fields}   
        for r in mi:
            if r['tsi1'] in sub_fields and r['tsi2'] in sub_fields2:
                mat[r['tsi1']][r['tsi2']] = float(r['mutual_info'])
    
        f, ax = plt.subplots(figsize=size)
        df = pandas.DataFrame(mat)
        sns.heatmap(
            df, square=True,vmin=0,vmax=1,
            xticklabels=fs.get_sparse_tsi(df.axes[1]),
            yticklabels=fs.get_sparse_tsi(df.axes[0])
        )
        plt.yticks(rotation='horizontal')
        plt.xticks(rotation='vertical')
        f.tight_layout()
        f.savefig(
            args.plot_path + 'tsi_heatmap_' + str(i) + '-' + str(j) + '.pdf'
        )
        plt.close()