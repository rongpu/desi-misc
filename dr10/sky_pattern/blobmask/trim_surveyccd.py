from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-decam-dr10-v2.fits'
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts', 'ra_bore', 'dec_bore', 'ra', 'dec', 'ccdraoff', 'ccddecoff', 'mjd_obs', 'plver']

ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey-ccds-decam-dr10-v2-trim.fits')

ccd = Table(fitsio.read(surveyccd_path))
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey-ccds-decam-dr10-v2-unique_exps.fits', overwrite=True)

##############################################################################################################

surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-decam-dr10-v2.fits'
ccd = Table(fitsio.read(surveyccd_path))
print(len(ccd))

# Remove zero-filled CCD entries
mask = ccd['expnum']!=0
ccd = ccd[mask]
print(len(ccd))

mask = ccd['filter']=='z'
ccd = ccd[mask]
print(len(ccd))

ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey-ccds-decam-dr10-v2-z_only.fits')