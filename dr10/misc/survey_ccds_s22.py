from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

tmp = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', columns=['ccdname']))
mask = tmp['ccdname']=='S22'
idx = np.where(mask)[0]
print(len(idx), len(idx)/len(tmp))

cat = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', rows=idx))
cat.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-s22.fits')

s22 = cat.copy()

# s22 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-s22.fits'))

##########################################################################################################################

cat = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits'))
_, idx = np.unique(cat['expnum'], return_index=True)
cat = cat[idx]
mask = ~np.in1d(cat['expnum'], s22['expnum'])
cat = cat[mask]
cat = vstack([s22, cat])

cat.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures.fits')
