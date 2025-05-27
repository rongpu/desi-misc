from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool


exp = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/survey-ccds-decam-dr11-trim.fits'))
print(len(exp))

mask = exp['filter']=='z'
exp = exp[mask]
print(len(exp))


def check_fringe_correction(index):
    img_fn = '/dvs_ro/cfs/cdirs/cosmo/staging/' + exp['image_filename'][index]
    if 'FRGSCALE' in fitsio.read_header(img_fn, ext=1).keys():
        old_fringe = True
    else:
        old_fringe = False
    if 'FRGSCNEW' in fitsio.read_header(img_fn, ext=1).keys():
        new_fringe = True
    else:
        new_fringe = False

    return old_fringe, new_fringe


n_processes = 256
with Pool(processes=n_processes) as pool:
    res = pool.map(check_fringe_correction, np.arange(len(exp)), chunksize=1)

old_fringe, new_fringe = np.array(res).T
exp['old_fringe'] = old_fringe
exp['new_fringe'] = new_fringe
exp.write('/global/cfs/cdirs/desi/users/rongpu/data/dr11dev/sky_pattern/survey-ccds-decam-dr11-unique_exps_z_fringe.fits')

old_fringe = exp['old_fringe']
new_fringe = exp['new_fringe']

mask = (old_fringe) & (~new_fringe)
exp1 = exp[mask].copy()
print(len(exp1))

with open("/global/cfs/cdirs/desi/users/rongpu/data/dr11dev/misc/fringe_list.txt", "w") as f:
    for fn in exp1['image_filename']:
        f.write(fn+'\n')
