# Compile a catalog of PSFEx profiles for a small number of exposures

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from multiprocessing import Pool

n_processes = 32
n_exposure_in_each_band = 200

psfex_columns = ['expnum', 'ccdname', 'psf_fwhm', 'psf_mask', 'moffat_alpha', 'moffat_beta', 'sum_diff', 'fit_original', 'failure']

output_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/calib/psfex'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'

ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'filter', 'image_filename', 'ccd_cuts']))
print(len(ccd))

mask = ccd['ccd_cuts']==0
ccd = ccd[mask]
print(len(ccd))

# Downsampling
np.random.seed(241)
expnum_list = []
for band in ['g', 'r', 'z']:
    mask = ccd['filter']==band
    expnum_list_all = np.unique(ccd['expnum'][mask])
    expnum_list += list(np.random.choice(expnum_list_all, size=n_exposure_in_each_band, replace=False))
print(len(expnum_list))


def get_psfex(expnum):

    mask = ccd['expnum']==expnum
    band = ccd['filter'][mask][0]

    image_filename = ccd['image_filename'][mask][0]
    psfex_filename = image_filename[:image_filename.find('.fits.fz')]+'-psfex.fits'
    psfex_path = os.path.join(output_dir, psfex_filename)

    if not os.path.isfile(psfex_path):
        raise ValueError('No PSFEx file found:', psfex_path)

    with fitsio.FITS(psfex_path) as hdu:
        if not 'moffat_alpha' in hdu[1].get_colnames():
            print('Error: No Moffat parameters found:', psfex_path, expnum)
            return None

    data = Table(fitsio.read(psfex_path, columns=psfex_columns))
    data['psf_mask0'] = data['psf_mask'][:, 0]
    data.remove_column('psf_mask')
    data['filter'] = band

    return data


def main():

    print('Starting')

    with Pool(processes=n_processes) as pool:
        res = pool.map(get_psfex, expnum_list)

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    psf_params = vstack(res)
    psf_params.write('/global/cfs/cdirs/desi/users/rongpu/dr9/psfex/survey-ccds-decam-dr9m-psfex.fits')

    print('Done!!!!!!!!!!!!!!!!!!!!!')


if __name__=="__main__":
    main()
