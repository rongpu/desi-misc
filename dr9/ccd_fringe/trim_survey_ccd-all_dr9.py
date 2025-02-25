# Trim survey-ccd table to save memory

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts', 'mjd_obs']
ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
# ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
# mask = ccd['ccd_cuts']==0
mask = ccd['filter']=='z' # include only z-band images
ccd = ccd[mask]

ccd.write('/global/project/projectdirs/desi/users/rongpu/dr9/fringe/temp/survey-ccds-decam-dr9-z-band-only-trim-all-dr9.fits')
