from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool


exp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey-ccds-decam-dr10-v2-unique_exps.fits'))
print(len(exp))

mask = exp['filter']=='z'
exp = exp[mask]
print(len(exp))


def check_fringe_correction(index):
    img_fn = '/global/project/projectdirs/cosmo/staging/' + exp['image_filename'][index]
    if 'FRGSCALE' in fitsio.read_header(img_fn, ext=1).keys():
        old_fringe = True
    else:
        old_fringe = False
    if 'FRGSCNEW' in fitsio.read_header(img_fn, ext=1).keys():
        new_fringe = True
    else:
        new_fringe = False

    return old_fringe, new_fringe


n_processes = 64
with Pool(processes=n_processes) as pool:
    res = pool.map(check_fringe_correction, np.arange(len(exp)), chunksize=1)

old_fringe, new_fringe = np.array(res).T
exp['old_fringe'] = old_fringe
exp['new_fringe'] = new_fringe
exp.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey-ccds-decam-dr10-v2-unique_exps_z_fringe.fits')
