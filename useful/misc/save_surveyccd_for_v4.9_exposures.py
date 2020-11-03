# Save the location of g-band images with ringing pupil artifacts

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

# survey-ccd table used in dr9m
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'
ccd = Table(fitsio.read(surveyccd_path, columns=['filter', 'ccdname', 'expnum', 'image_filename', 'ccd_cuts', 'plver']))
print(len(ccd))
# Only keep unique exposures
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
print(len(ccd))

# the older version before the V4.9 replacement
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
ccd1 = Table(fitsio.read(surveyccd_path, columns=['filter', 'ccdname', 'expnum', 'image_filename', 'ccd_cuts', 'plver']))
print(len(ccd1))
# Only keep unique exposures
_, idx = np.unique(ccd1['expnum'], return_index=True)
ccd1 = ccd1[idx]
print(len(ccd1))

mask = ccd['plver']=='V4.9'
print(np.sum(mask))
expnum_list = ccd['expnum'][mask]
mask1 = np.in1d(ccd1['expnum'], expnum_list)

ccd = ccd[mask]
ccd1 = ccd1[mask1]

ccd.write('/global/cfs/cdirs/desi/users/rongpu/tmp/survey-ccds-decam-dr9m.fits')
ccd1.write('/global/cfs/cdirs/desi/users/rongpu/tmp/survey-ccds-decam-dr9.fits')
