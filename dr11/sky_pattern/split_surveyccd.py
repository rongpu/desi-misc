from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr11-early/survey-ccds-decam-dr11-merged-incl-early.fits'

ccd = Table(fitsio.read(surveyccd_path))
print(len(ccd))

# only do z-band sky templates
mask = ccd['filter']=='z'
ccd = ccd[mask]
print(len(ccd))

np.random.seed(393)
expnum_list_all = np.unique(ccd['expnum'])
expnum_list_all = np.random.choice(expnum_list_all, size=len(expnum_list_all), replace=False)

n_chunks = 20
expnum_list_split = np.array_split(expnum_list_all, n_chunks)

for index, expnum_list in enumerate(expnum_list_split):
    mask = np.in1d(ccd['expnum'], expnum_list)
    ccd1 = ccd[mask].copy()
    ccd1.write('/global/cfs/cdirs/desi/users/rongpu/data/dr11dev/sky_pattern/survey_ccds_chunks/survey-ccds-decam-dr11-merged-incl-early-{}.fits'.format(index))

# Sanity check
ccd_id_str = []
for index in range(n_chunks):
    tmp = fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/dr11dev/sky_pattern/survey_ccds_chunks/survey-ccds-decam-dr11-merged-incl-early-{}.fits'.format(index), columns=['expnum', 'ccdname'])
    ccd_id_str.append(np.char.add(np.array(tmp['expnum']).astype(str), tmp['ccdname']))
ccd_id_str = np.concatenate(ccd_id_str)
ccd_id_str.sort()
ccd['ccd_id_str'] = np.char.add(np.array(ccd['expnum']).astype(str), ccd['ccdname'])
ccd.sort('ccd_id_str')
assert len(ccd_id_str)==len(ccd) and np.all(ccd_id_str==ccd['ccd_id_str'])
