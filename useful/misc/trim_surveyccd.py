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
