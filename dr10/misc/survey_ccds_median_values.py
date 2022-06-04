from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool


n_processes = 128

columns = ['expnum', 'filter', 'ccdname', 'fwhm', 'skyrms', 'sig1', 'ccdskycounts', 'ccdskysb', 'ccd_cuts']
median_columns = ['fwhm', 'skyrms', 'sig1', 'ccdskycounts', 'ccdskysb']
median_columns1 = ['median_fwhm', 'median_skyrms', 'median_sig1', 'median_ccdskycounts', 'median_ccdskysb']

cat = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', columns=columns))
print(len(cat))


def get_medians(exp_idx):

    expnum_list = expnum_unique[exp_idx]

    tt = Table()
    tt['expnum'] = expnum_list

    arr = np.zeros([len(expnum_list), len(median_columns1)])
    tt = hstack([tt, Table(arr, names=median_columns1)])
    tt['n_ccds'] = 0
    tt['n_good_ccds'] = 0

    for index in np.arange(len(exp_idx)):

        idx = exp_order[exp_counts[exp_idx[index]]:exp_counts[exp_idx[index]+1]]

        for ii in range(len(median_columns)):
            colname_old = median_columns[ii]
            colname_new = median_columns1[ii]
            tt[colname_new][index] = np.nanmedian(cat[colname_old][idx])

        tt['n_ccds'][index] = len(idx)
        band = cat['filter'][idx][0]
        if band!='Y':
            mask = (cat['ccd_cuts'][idx]==0) | (cat['ccd_cuts'][idx]==2**14)
            tt['n_good_ccds'][index] = np.sum(mask)
        else:
            ccd_cuts = cat['ccd_cuts'][idx].copy()
            ccd_cuts[ccd_cuts&2**1>0] -= 2**1  # ignore NOT_GRZ
            ccd_cuts[ccd_cuts&2**15>0] -= 2**15  # ignore TOO_MANY_BAD_CCDS
            mask = (ccd_cuts==0) | (ccd_cuts==2**14)
            tt['n_good_ccds'][index] = np.sum(mask)
    return tt


# pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)
expnum_unique, exp_counts = np.unique(cat['expnum'], return_counts=True)

exp_counts = np.insert(exp_counts, 0, 0)
exp_counts = np.cumsum(exp_counts)

exp_order = np.argsort(cat['expnum'])

# split among the Cori processors
exp_idx_split = np.array_split(np.arange(len(expnum_unique)), n_processes)

# start multiple worker processes
with Pool(processes=n_processes) as pool:
    res = pool.map(get_medians, exp_idx_split)

cat_stack = vstack(res)
cat_stack.sort('expnum')

cat_stack.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures-medians.fits', overwrite=True)

