# Add DR9 z-band exposures that now require sky fitting

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

ccd = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/survey-ccds-decam-dr11-trim.fits'))
print('CCDs', len(ccd))

mask = ccd['filter']=='z'
ccd = ccd[mask]
print('CCDs', len(ccd))

dr9_reuse = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/sky_pattern/sky-scales-reuse_dr9.fits'))
skyrun_dr11 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/skyrunsdr11-v3-rearrange.fits'))

mask = (~np.in1d(ccd['expnum'], dr9_reuse['expnum'])) & (ccd['mjd_obs']<skyrun_dr11['mjd_obs'].min())
ccd = ccd[mask]
print('CCDs', len(ccd))

# from notebooks/dr11/exposures_with_high_stellar_density.ipynb
fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/high_stellar_density_expnum_list.txt'
with open(fn, 'r') as f:
    expnum_list = f.read().splitlines()
expnum_list = [int(ii) for ii in expnum_list]
print(len(expnum_list))
mask = (~np.in1d(ccd['expnum'], expnum_list))
ccd = ccd[mask]
print('CCDs', len(ccd))

######################### Get DR9 template indices #########################
dr9 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr9/sky_pattern/skyrunsgoodcountexpnumv48dr8.fits'))
print('dr9', len(dr9))

mask = dr9['filter']=='z'
dr9 = dr9[mask]
print('dr9', len(dr9))

ccd['run'] = -1
for index in range(len(ccd)):
    ii = np.argmin(np.abs(dr9['mjd_obs']-ccd['mjd_obs'][index]))
    ccd['run'][index] = dr9['run'][ii]

ccd = ccd[['expnum', 'run', 'filter', 'mjd_obs', 'image_filename']]
skyrun_dr11 = skyrun_dr11[['expnum', 'run', 'filter', 'mjd_obs', 'image_filename']]
skyrun_dr11 = vstack([ccd, skyrun_dr11], join_type='exact')

skyrun_dr11.sort('expnum')

skyrun_dr11.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/skyrunsdr11-v3-rearrange-dr9_added.fits', overwrite=False)
