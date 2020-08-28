from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from multiprocessing import Pool

n_processes = 32

output_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/calib/psfex'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'

ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'filter', 'image_filename']))

expnum_list = np.unique(ccd['expnum'])
print(len(expnum_list))

# # Downsampling
# expnum_list = np.sort(np.random.choice(expnum_list, size=64, replace=False))
# print(len(expnum_list))

def get_moffat_params(expnum):

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

        hdu = fits.open(psfex_path)
        data = Table(hdu[1].data)[['expnum', 'ccdname', 'plver', 'psf_patch_ver', 'moffat_alpha', 'moffat_beta', 'sum_diff', 'fit_original', 'failure']]
        data['filter'] = band

        return data

def main():

    print('Starting')

    with Pool(processes=n_processes) as pool:
        res = pool.map(get_moffat_params, expnum_list)

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    psf_params = vstack(res)
    psf_params.write('/global/homes/r/rongpu/data/survey-ccds-decam-dr9m-psfex-moffat-params.fits')


    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

