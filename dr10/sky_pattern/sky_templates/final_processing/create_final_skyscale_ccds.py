from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits


ss = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_scales/skyscales_ccds_raw.fits'))
print(len(ss))

t = Table()
t['expnum'], t['count'] = np.unique(ss['expnum'], return_counts=True)
t.sort('count')
# t

# Require 50 CCDs with valid skyscale fits
mask = t['count']>=50
np.sum(mask)/len(mask)
good_exp = t['expnum'][mask]
mask = np.in1d(ss['expnum'], good_exp)
ss = ss[mask]

_, idx = np.unique(ss['expnum'], return_index=True)
ss = ss[idx]
print(len(ss))

ss.remove_columns(['image_hdu', 'ccdname', 'ccdskyscale'])
ss.rename_column('medianskyscale', 'skyscale')

# add plprocid
ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey-ccds-decam-dr10-v2-unique_exps.fits', columns=['expnum', 'plprocid']))
ss = join(ss, ccd, keys='expnum')
print(len(ss))

ss.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_scales/skyscales.fits', overwrite=True)

