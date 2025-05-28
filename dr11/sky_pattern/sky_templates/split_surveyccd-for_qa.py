from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits


# skyrun = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/skyrunsdr11-v3-rearrange-dr9_added.fits'))
# skyrun = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/schlafly/decals/skyrunsdr11-v3.fits'))
skyrun = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/dr11dev/sky_scales/skyscales.fits'))
print('skyrun', len(skyrun))

# mask = skyrun['filter']=='z'
# skyrun = skyrun[mask]
# print('skyrun', len(skyrun))

surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr11-early/survey-ccds-decam-dr11-merged-incl-early.fits'
ccd = Table(fitsio.read(surveyccd_path))
print('ccd', len(ccd))

mask = np.in1d(ccd['expnum'], skyrun['expnum'])
ccd = ccd[mask]
print('ccd', len(ccd))

run_list = np.unique(skyrun['run'])

for run in run_list:
    mask = skyrun['run']==run
    expnum_list = np.array(skyrun['expnum'][mask])
    mask = np.in1d(ccd['expnum'], expnum_list)
    ccd1 = ccd[mask].copy()
    fn = '/pscratch/sd/r/rongpu/tmp/survey_ccds_in_runs_for_qa/survey-ccds-decam-dr11-merged-incl-early-run_{}.fits'.format(run)
    ccd1.write(fn)

# Sanity check
ccd_id = []
for run in run_list:
    tmp = fitsio.read('/pscratch/sd/r/rongpu/tmp/survey_ccds_in_runs_for_qa/survey-ccds-decam-dr11-merged-incl-early-run_{}.fits'.format(run), columns=['expnum', 'image_hdu'])
    ccd_id.append(np.array(100*tmp['expnum'] + tmp['image_hdu']))
ccd_id = np.concatenate(ccd_id)
ccd_id.sort()
ccd['ccd_id'] = 100*ccd['expnum'] + ccd['image_hdu']
ccd.sort('ccd_id')
assert len(ccd_id)==len(ccd) and np.all(ccd_id==ccd['ccd_id'])
    