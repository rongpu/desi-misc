from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

from multiprocessing import Pool
from scipy import ndimage


surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr11-early/survey-ccds-decam-dr11.fits'
# surveyccd_path = '/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/dr11/survey-ccds-decam-dr11-merged.fits'
# surveyccd_path = '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr10/survey-ccds-decam-dr10.fits.gz'
# surveyccd_path = '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'

ccd = Table(fitsio.read(surveyccd_path, columns=['ccdname', 'expnum']))
print(len(ccd))

ccd['idx'] = np.arange(len(ccd))

# Unique exposures
print(len(ccd), len(np.unique(ccd['expnum'])))
# prefer N1, N7, S1 (all good edge CCDs):
ccd['priority'] = 100
for ii, ccdname in enumerate(['N1', 'N7', 'S1']):
    mask = ccd['ccdname']==ccdname
    ccd['priority'][mask] = ii
# deprioritize potentially bad CCDs
ccd_blacklist = ['N10', 'N13', 'N15', 'N30', 'S30', 'S7 ', 'S7']
mask = np.in1d(ccd['ccdname'], ccd_blacklist)
ccd['priority'][mask] = 200

ccd.sort('priority')
_, idx_keep = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx_keep]
ccd.remove_column('priority')
print(len(ccd), len(np.unique(ccd['expnum'])))

idx = np.sort(ccd['idx'])
ccd = Table(fitsio.read(surveyccd_path, rows=idx))

ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11/survey-ccds-decam-dr11-merged-trim-new.fits')
# ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11/survey-ccds-decam-dr11-merged-trim.fits')
# ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11/survey-ccds-decam-dr10-trim.fits')
# ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11/survey-ccds-decam-dr9-trim.fits')
