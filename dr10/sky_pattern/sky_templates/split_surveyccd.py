from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits


skyrun = Table.read('/global/cfs/cdirs/desi/users/schlafly/decals/skyrunsdr10-v3.fits')

surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-decam-dr10-v2.fits'
ccd = Table(fitsio.read(surveyccd_path))
print(len(ccd))

# Remove zero-filled CCD entries
mask = ccd['expnum']!=0
ccd = ccd[mask]
print(len(ccd))

mask = np.in1d(ccd['expnum'], skyrun['expnum'])
ccd = ccd[mask]
print(len(ccd))

run_list = np.unique(skyrun['run'])

for run in run_list:
    mask = skyrun['run']==run
    expnum_list = np.array(skyrun['expnum'][mask])
    mask = np.in1d(ccd['expnum'], expnum_list)
    ccd1 = ccd[mask].copy()
    ccd1.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey_ccds_in_runs/survey-ccds-decam-dr10-v2-run_{}.fits'.format(run))

# Sanity check
ccd_id = []
for run in run_list:
    tmp = fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey_ccds_in_runs/survey-ccds-decam-dr10-v2-run_{}.fits'.format(run), columns=['expnum', 'image_hdu'])
    ccd_id.append(np.array(100*tmp['expnum'] + tmp['image_hdu']))
ccd_id = np.concatenate(ccd_id)
ccd_id.sort()
ccd['ccd_id'] = 100*ccd['expnum'] + ccd['image_hdu']
ccd.sort('ccd_id')
print(len(ccd_id)==len(ccd), np.all(ccd_id==ccd['ccd_id']))
