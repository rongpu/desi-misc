from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack

columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'mjd_obs', 'ra', 'dec', 'ccdraoff', 'ccddecoff', 'ccd_cuts', 'plver', 'ccdskycounts']

ccd = Table.read('/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz')
ccd = ccd[columns]

for band in ['g', 'r', 'z']:
    mask = ccd['filter']==band
    ccd1 = (ccd[mask]).copy()
    ccd1.write('/global/cfs/cdirs/desi/users/rongpu/dr9/misc/survey-ccds-decam-dr9-trim-{}.fits'.format(band))
