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

# surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
# ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'image_filename']))
# print(len(ccd))
# _, idx = np.unique(ccd['expnum'], return_index=True)
# ccd = ccd[idx]
# print(len(ccd))
# ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-old-trim.fits')

# surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'
# ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'image_filename']))
# print(len(ccd))
# _, idx = np.unique(ccd['expnum'], return_index=True)
# ccd = ccd[idx]
# print(len(ccd))
# ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-trim.fits')

surveyccd_path = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'image_filename']))
print(len(ccd))
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
print(len(ccd))
ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-final-trim.fits')

surveyccd_path = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'image_filename', 'ccdname', 'image_hdu', 'filter']))
print(len(ccd))
mask = ccd['filter']=='z'
ccd = ccd[mask]
print(len(ccd))
ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-final-trim-z.fits')

surveyccd_path_dr8 = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9-garage/reorg/decam/survey-ccds-decam-dr8-newlocs2.fits.gz'
ccd = Table(fitsio.read(surveyccd_path_dr8, columns=['expnum', 'image_filename']))
print(len(ccd))
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
print(len(ccd))
ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr8-trim.fits')

#####################################

camera_list = ['mosaic', '90prime', 'decam']
for camera in camera_list:
    surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-{}-dr9.fits.gz'.format(camera)
    columns = ['image_filename', 'image_hdu', 'expnum', 'plver', 'plprocid', 'ccdname', 'propid', 'filter', 'exptime', 'mjd_obs', 'ra_bore', 'dec_bore', 'zpt', 'ccdskycounts', 'ccdskysb', 'ccdphrms', 'ccdnphotom', 'ccd_cuts']
    ccd = Table(fitsio.read(surveyccd_path, columns=columns))
    print(len(ccd))
    # Add ccd_cuts_ok and median values
    ccd['ccd_cuts_ok'] = False
    ccd['median_ccdskycounts'] = 0.
    ccd['median_ccdskysb'] = 0.
    ccd['median_ccdphrms'] = 0.
    expnum_list = np.unique(ccd['expnum'])
    for ii, expnum in enumerate(expnum_list):
        if ii%1000==0:
            print('{}/{}'.format(ii, len(expnum_list)))
        mask = ccd['expnum']==expnum
        ccd['median_ccdskycounts'][mask] = np.median(ccd['ccdskycounts'][mask])
        ccd['median_ccdskysb'][mask] = np.median(ccd['ccdskysb'][mask])
        ccd['median_ccdphrms'][mask] = np.median(ccd['ccdphrms'][mask])
        if (np.sum(ccd['ccd_cuts'][mask]==0)) >= (np.sum(mask)//2):
            ccd['ccd_cuts_ok'][mask] = True
    _, idx = np.unique(ccd['expnum'], return_index=True)
    ccd = ccd[idx]
    print(len(ccd))
    ccd.write('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-{}-dr9-trim-less.fits'.format(camera), overwrite=True)
