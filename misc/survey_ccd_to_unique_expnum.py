from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits


# columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts', 'ra', 'dec', 'ra_bore', 'dec_bore', 'plver']

ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'))
print(len(ccd))
_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx]
print(len(exp))
exp.write('/global/cscratch1/sd/rongpu/temp/survey-ccds-decam-dr9-unique-exp.fits')


ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-decam-dr10.fits'))
print(len(ccd))
_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx]
print(len(exp))
exp.write('/global/cscratch1/sd/rongpu/temp/survey-ccds-decam-dr10-20211119-unique-exp.fits')

