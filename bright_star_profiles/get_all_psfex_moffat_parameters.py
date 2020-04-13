from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

output_dir = '/global/cfs/projectdirs/cosmo/work/legacysurvey/dr9/calib/patched-psfex'
surveyccd_path = '/global/cfs/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'

ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
psf_params_stack = []

unique_expnum = np.unique(ccd['expnum'])
print(len(unique_expnum))

# Downsampling
unique_expnum = np.sort(np.random.choice(unique_expnum, size=2000, replace=False))

for index, expnum in enumerate(unique_expnum):

    if index%100==0:
        print(index, '/', len(unique_expnum))

    mask = ccd['expnum']==expnum
    band = ccd['filter'][mask][0]

    image_filename = ccd['image_filename'][mask][0]
    psfex_filename_new = image_filename[:image_filename.find('.fits.fz')]+'-psfex.fits'
    psfex_path_new = os.path.join(output_dir, psfex_filename_new)

    if not os.path.isfile(psfex_path_new):
        raise ValueError

    with fitsio.FITS(psfex_path_new) as hdu:
        if not 'moffat_alpha' in hdu[1].get_colnames():
            continue

        hdu = fits.open(psfex_path_new)
        data = Table(hdu[1].data)[['expnum', 'ccdname', 'psf_patch_ver', 'moffat_alpha', 'moffat_beta', 'sum_diff', 'fit_original', 'failure']]
        data['filter'] = band

        psf_params_stack.append(data)

psf_params_stack = vstack(psf_params_stack)
psf_params_stack.write('/global/homes/r/rongpu/data/survey-ccds-decam-dr9-cut-psfex-moffat-params-subsample.fits')
