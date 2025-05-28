from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits


ss = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/sky_scales/skyscales_ccds_raw.fits'))
print(len(ss))

ss.remove_columns(['image_hdu', 'ccdname', 'ccdskyscale'])
ss.rename_column('medianskyscale', 'skyscale')

# Require 50 CCDs with valid skyscale fits
t = Table()
t['expnum'], t['count'] = np.unique(ss['expnum'], return_counts=True)
mask = t['count']>=50
print(np.sum(mask)/len(mask))
good_exp = t['expnum'][mask]
mask = np.in1d(ss['expnum'], good_exp)
ss = ss[mask]

_, idx = np.unique(ss['expnum'], return_index=True)
ss = ss[idx]
print(len(ss))

# add plprocid
expnum_arr = np.array(ss['expnum'])
ccd = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/survey-ccds-decam-dr11-trim.fits', columns=['expnum', 'plprocid', 'filter']))
ss = join(ss, ccd[['expnum', 'plprocid']], keys='expnum', join_type='left')
print(len(ss))
assert np.all(ss['expnum']==expnum_arr)

ss_dr9 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/sky_scales/sky-scales-reuse_dr9.fits'))
ss_dr9['dr9_skyscale'] = True
ss['dr9_skyscale'] = False
ss = vstack([ss_dr9, ss])
print(len(ss))

# add filter
ss = join(ss, ccd[['expnum', 'filter']], keys='expnum', join_type='left')
print(len(ss))

ss.sort('expnum')
ss.sort('run', kind='stable')

ss.write('/global/cfs/cdirs/desi/users/rongpu/data/dr11dev/sky_scales/skyscales.fits', overwrite=False)

