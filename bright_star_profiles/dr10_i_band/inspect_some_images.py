# srun -N 1 -C haswell -c 64 -t 04:00:00 -q interactive python inspect_some_images.py

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

from pathlib import Path
from multiprocessing import Pool

from decam_postage_stamps import decam_plot


def make_a_plot(img_path):

    plot_path = os.path.basename(img_path).replace('.fits.fz', '.png')
    plot_path = os.path.join('plots/DR10b', plot_path)
    print(plot_path)

    decam_plot(img_path, plot_path, blob_mask=False)


time_start = time.time()

n_processes = 32

img_list = glob.glob('/global/cfs/cdirs/cosmo/staging/decam/DECam_CP-DR10b/*/*_ooi_*.fits.fz')
print(len(img_list))

np.random.seed(8721)
img_list = np.random.choice(img_list, 128, replace=False)

with Pool(processes=n_processes) as pool:
    pool.map(make_a_plot, img_list)

print('All done!', time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
