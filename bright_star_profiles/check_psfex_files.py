from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

output_dir = '/global/cfs/projectdirs/cosmo/work/legacysurvey/dr9/calib/patched-psfex'
# output_dir = '/global/cfs/projectdirs/cosmo/work/legacysurvey/dr9/calib/v2-patched-psfex'
surveyccd_path = '/global/cfs/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'

ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)

ccds_to_reprocess = []

unique_expnum = np.unique(ccd['expnum'])
print(len(unique_expnum))

for index, expnum in enumerate(unique_expnum):

    if index%100==0:
        print(index, '/', len(unique_expnum))

    mask = ccd['expnum']==expnum
    band = ccd['filter'][mask][0]
    ccd_index = np.where(mask)[0][0]

    image_filename = ccd['image_filename'][mask][0]
    psfex_filename_new = image_filename[:image_filename.find('.fits.fz')]+'-psfex.fits'
    psfex_path_new = os.path.join(output_dir, psfex_filename_new)

    if not os.path.isfile(psfex_path_new):
        raise ValueError

    with fitsio.FITS(psfex_path_new) as hdu:
        if not 'moffat_alpha' in hdu[1].get_colnames():
            ccds_to_reprocess.append(ccd_index)
            # print('found one: expnum =', expnum)
            print(ccd['image_filename'][ccd_index])

if len(ccds_to_reprocess)>0:
    ccd = ccd[ccds_to_reprocess]
    # ccd.write('/global/homes/r/rongpu/data/survey-ccds-decam-dr9-cut-no_moffat_params-v2.fits')
    ccd.write('/global/homes/r/rongpu/data/survey-ccds-decam-dr9-cut-no_moffat_params.fits')
