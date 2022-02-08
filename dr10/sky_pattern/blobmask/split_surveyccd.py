from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-decam-dr10-v2.fits'

ccd = Table(fitsio.read(surveyccd_path))
print(len(ccd))

# Remove zero-filled CCD entries
mask = ccd['expnum']!=0
ccd = ccd[mask]
print(len(ccd))

np.random.seed(393)
expnum_list_all = np.unique(ccd['expnum'])
expnum_list_all = np.random.choice(expnum_list_all, size=len(expnum_list_all), replace=False)

n_chunks = 200
expnum_list_split = np.array_split(expnum_list_all, n_chunks)

for index, expnum_list in enumerate(expnum_list_split):
    mask = np.in1d(ccd['expnum'], expnum_list)
    ccd1 = ccd[mask].copy()
    ccd1.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey_ccds_chunks/survey-ccds-decam-dr10-v2-{}.fits'.format(index))


# Sanity check
ccd_id = []
for index in range(n_chunks):
    tmp = fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey_ccds_chunks/survey-ccds-decam-dr10-v2-{}.fits'.format(index), columns=['expnum', 'image_hdu'])
    ccd_id.append(np.array(100*tmp['expnum'] + tmp['image_hdu']))
ccd_id = np.concatenate(ccd_id)
ccd_id.sort()
ccd['ccd_id'] = 100*ccd['expnum'] + ccd['image_hdu']
ccd.sort('ccd_id')
print(len(ccd_id)==len(ccd), np.all(ccd_id==ccd['ccd_id']))
