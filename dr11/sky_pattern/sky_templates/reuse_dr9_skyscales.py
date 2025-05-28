# Reuse DR9 skyscales for some exposures

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

skyscales = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/calib/sky_pattern/sky-scales.fits'))
print('skyscales', len(skyscales))

_, idx = np.unique(skyscales['expnum'], return_index=True)
skyscales = skyscales[idx]
print('skyscales', len(skyscales))

skyscales.sort('expnum')
skyscales.sort('run')

skyscales.sort('expnum')
skyscales.sort('run', kind='stable')

dr11 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/survey-ccds-decam-dr11-trim.fits'))
print('dr11', len(dr11))

dr9 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/survey-ccds-decam-dr9-trim.fits'))
print('dr9', len(dr9))

idx = np.intersect1d(dr11['expnum'], dr9['expnum'])
mask = np.in1d(dr11['expnum'], idx)
dr11 = dr11[mask]
mask = np.in1d(dr9['expnum'], idx)
dr9 = dr9[mask]
dr11.sort('expnum')
dr9.sort('expnum')
assert np.all(dr11['expnum']==dr9['expnum'])
print('dr11', len(dr11))
print('dr9', len(dr9))

# edge glow
run_list = [446, 436, 437, 438, 439, 445, 447, 453, 455, 456, 469, 471]

mask = np.in1d(skyscales['run'], run_list)
mask |= skyscales['filter']=='z'
skyscales = skyscales[mask]
print('skyscales', len(skyscales))

# Only reuse skyscales for identical CP versions
mask1 = (dr11['plver']==dr9['plver']) & (dr11['plprocid']==dr9['plprocid'])
mask = np.in1d(skyscales['expnum'], dr11['expnum'][mask1])
skyscales = skyscales[mask]
print('skyscales', len(skyscales))

skyscales.rename_column('PLPROCID', 'plprocid')
skyscales = skyscales[['skyscale', 'expnum', 'run', 'plprocid']]

skyscales.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/sky_scales/sky-scales-reuse_dr9.fits')
