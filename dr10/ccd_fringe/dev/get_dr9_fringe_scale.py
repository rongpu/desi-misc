from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool


exp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/misc/survey-ccds-dr10-deep-fields-v1-z-fringe-stacking.fits'))
print(len(exp))

_, idx = np.unique(exp['expnum'], return_index=True)
exp = exp[idx]
print(len(exp))

mask = exp['filter']=='z'
exp = exp[mask]
print(len(exp))


def get_frgscale(index):
    img_fn = '/global/cfs/cdirs/cosmo/staging/' + exp['image_filename'][index]
    nccd = len(fitsio.FITS(img_fn))-1

    cat = Table()
    cat['frgscale_old'] = np.zeros(nccd)
    cat['frgscale_new'] = np.zeros(nccd)
    cat['ccdname'] = '   '

    for ii, hdu_index in enumerate(range(1, nccd+1)):
        header = fitsio.read_header(img_fn, ext=hdu_index)
        cat['ccdname'][ii] = header['DETPOS']
        if 'FRGSCALE' in header.keys():
            cat['frgscale_old'][ii] = header['FRGSCALE']
        else:
            cat['frgscale_old'][ii] = -99.
        if 'FRGSCNEW' in header.keys():
            cat['frgscale_new'][ii] = header['FRGSCNEW']
        else:
            cat['frgscale_new'][ii] = -99.

    mask = cat['frgscale_new']!=-99
    cat['n_fringe_new'] = np.sum(mask)
    cat['frgscale_new_median'] = np.median(cat['frgscale_new'][mask])
    cat['expnum'] = exp['expnum'][index]

    return cat


n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(get_frgscale, np.arange(len(exp)), chunksize=1)

cat = vstack(res)

cat.write('/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-dr9-frgscale.fits', overwrite=True)
