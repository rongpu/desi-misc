from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/useful/'))
from decam_postage_stamps import decam_plot


n_process = 128
image_dir = '/global/cfs/cdirs/cosmo/staging'


def wrapper(index):
    decam_plot(os.path.join(image_dir, exp['image_filename'][index]),
               plot_path='/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_postage/z-fringe-stack/{}/{}.png'.format(exp['plver'][index], exp['expnum'][index]))


cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/misc/survey-ccds-dr10-deep-fields-v1-z-fringe-stacking.fits'))
_, idx = np.unique(cat['expnum'], return_index=True)
exp = cat[idx]

mask_keep = np.full(len(exp), False)
for plver in np.unique(exp['plver']):
    mask = exp['plver']==plver
    idx = np.where(mask)[0]
    if np.sum(mask)>20:
        np.random.seed(38123)
        idx = np.random.choice(idx, size=20, replace=False)
    mask_keep[idx] = True
exp = exp[mask_keep]
print(len(exp))

with Pool(processes=n_process) as pool:
    res = pool.map(wrapper, np.arange(len(exp)))
