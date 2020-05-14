from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits

skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8.fits')

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
ccd_columns = ['image_hdu', 'expnum', 'ccdname', 'filter', 'ccdskycounts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))

ccd['ccd_id'] = 100*ccd['expnum'] + ccd['image_hdu']

skyscale = Table.read('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_raw.fits')
skyscale.remove_column('ccdname')

ccd_new = join(ccd, skyscale, keys=['image_hdu', 'expnum'], join_type='left')
print(len(ccd))
print(len(ccd_new))

mask = np.copy(ccd_new['ccdskyscale'].mask)
ccd_new = ccd_new.filled(fill_value=0)
ccd_new['ccdskyscale'][mask] = 0
ccd_new['medianskyscale'][mask] = 0
ccd_new['run'][mask] = -1

# Line match to survey-ccds
ccd_reverse_sort = np.array(ccd['ccd_id']).argsort().argsort()
ccd_new.sort('ccd_id')
ccd_new = ccd_new[ccd_reverse_sort]

ccd_new['skyscale'] = 0.

for index in range(len(skyrun)):

    if index%100==0:
        print(index, len(skyrun))

    expnum = skyrun['expnum'][index]
    run = skyrun['run'][index]

    # All CCDs in the exposure
    mask_exp = (ccd_new['expnum']==expnum)

    # CCDs with valid skyscales measured
    mask_has_skyscale = (ccd_new['run'][mask_exp]>=0)

    if np.sum(mask_has_skyscale)>=10:
        ccd_new['skyscale'][mask_exp] = (ccd_new['medianskyscale'][mask_exp])[mask_has_skyscale][0]
    else:
        ccd_new['skyscale'][mask_exp] = np.median(ccd_new['ccdskycounts'][mask_exp])

    # Fill in the run numbers last
    ccd_new['run'] = run

ccd_new.write('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_debug.fits')

# Remove the skyscale values that won't be used
ccd_new.remove_columns(['ccdskyscale', 'medianskyscale', 'ccdskycounts', 'ccd_id'])

ccd_new.write('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds.fits')
