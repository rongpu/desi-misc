# Create trimmed surveyccd table for decam_focal_plane_plot.py

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'image_filename']))
print(len(ccd))
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
print(len(ccd))
ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-old-trim.fits')

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'
ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'image_filename']))
print(len(ccd))
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
print(len(ccd))
ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-trim.fits')

surveyccd_path_dr8 = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9-garage/reorg/decam/survey-ccds-decam-dr8-newlocs2.fits.gz'
ccd = Table(fitsio.read(surveyccd_path_dr8, columns=['expnum', 'image_filename']))
print(len(ccd))
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
print(len(ccd))
ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr8-trim.fits')

#####################################

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'
columns = ['image_filename', 'image_hdu', 'expnum', 'plver', 'plprocid', 'ccdname', 'propid', 'filter', 'exptime', 'mjd_obs', 'ra_bore', 'dec_bore', 'ccd_cuts', 'ccdzpt']
ccd = Table(fitsio.read(surveyccd_path, columns=columns))
print(len(ccd))
# Add ccd_cuts_ok and median_ccdzpt
ccd['ccd_cuts_ok'] = False
ccd['median_ccdzpt'] = 0.
expnum_list = np.unique(ccd['expnum'])
for ii, expnum in enumerate(expnum_list):
    if ii%1000==0:
        print('{}/{}'.format(ii, len(expnum_list)))
    mask = ccd['expnum']==expnum
    if (np.sum(ccd['ccd_cuts'][mask]==0)) >= (np.sum(mask)//2):
        ccd['ccd_cuts_ok'][mask] = True
    ccd['median_ccdzpt'][mask] = np.median(ccd['ccdzpt'][mask])
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
print(len(ccd))
ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-trim-less.fits', overwrite=True)
