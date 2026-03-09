# srun -N 1 -C cpu -c 256 -t 04:00:00 -q interactive python get_all_psfex_moffat_parameters.py > get_all_psfex_moffat_parameters.log

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits

from multiprocessing import Pool


def get_moffat_params(expnum):

    mask = ccd['expnum']==expnum
    band = ccd['filter'][mask][0]

    image_filename = ccd['image_filename'][mask][0]
    psfex_filename = image_filename[:image_filename.find('.fits.fz')]+'-psfex.fits'
    psfex_path = os.path.join(psfex_dir, psfex_filename)

    if not os.path.isfile(psfex_path):
        raise ValueError('No PSFEx file found:', psfex_path)

    # with fitsio.FITS(psfex_path) as hdu:
    #     if not 'moffat_alpha' in hdu[1].get_colnames():
    #         print('Error: No Moffat parameters found:', psfex_path, expnum)
    #         return None

    # data = Table(hdu[1].data)[['expnum', 'ccdname', 'plver', 'psf_patch_ver', 'moffat_alpha', 'moffat_beta', 'sum_diff', 'fit_original', 'failure']]
    data = Table(fitsio.read(psfex_path))
    data['filter'] = band
    if 'psf_patch_ver' in data.colnames:
        data = data[['expnum', 'ccdname', 'filter', 'plver', 'psf_fwhm', 'psf_patch_ver', 'moffat_alpha', 'moffat_beta', 'sum_diff', 'fit_original', 'failure']]
    else:
        data = data[['expnum', 'ccdname', 'filter', 'plver', 'psf_fwhm']]

    return data


time_start = time.time()
print('Starting')

n_processes = 256

psfex_dir = '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/calib/psfex'
surveyccd_path = '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'

ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'filter', 'image_filename', 'ccdname', 'plver', 'ccd_cuts']))
print(len(ccd))

ccd['ccd_id_str'] = np.char.add(np.array(ccd['expnum']).astype(str), ccd['ccdname'])

expnum_list = np.unique(ccd['expnum'])
print(len(expnum_list))

# # Downsampling
# expnum_list = np.sort(np.random.choice(expnum_list, size=64, replace=False))
# print(len(expnum_list))

with Pool(processes=n_processes) as pool:
    res = pool.map(get_moffat_params, expnum_list, chunksize=1)

# # Remove None elements from the list
# for index in range(len(res)-1, -1, -1):
#     if res[index] is None:
#         res.pop(index)

psf_params = vstack(res).filled(-99)
print(len(psf_params))

psf_params['ccd_id_str'] = np.char.add(np.array(psf_params['expnum']).astype(str), psf_params['ccdname'])

# remove CCDs not in the CCD list
mask = np.in1d(psf_params['ccd_id_str'], ccd['ccd_id_str'])
psf_params = psf_params[mask]
print(len(psf_params))

# -99 fill CCDs missing from the CCD list
mask = ~np.in1d(ccd['ccd_id_str'], psf_params['ccd_id_str'])
psf_params = vstack([psf_params, ccd[['ccd_id_str', 'expnum', 'ccdname', 'filter', 'plver']][mask]], join_type='outer').filled(-99)
print(len(psf_params))

# Matching psf_params to ccd
if len(ccd)!=len(psf_params) or not np.all(np.unique(ccd['ccd_id_str'])==np.unique(psf_params['ccd_id_str'])):
    raise ValueError('ccd and psf_params have different id list')
t1_reverse_sort = np.array(ccd['ccd_id_str']).argsort().argsort()
psf_params = psf_params[np.argsort(psf_params['ccd_id_str'])[t1_reverse_sort]]
assert np.all(ccd['ccd_id_str']==psf_params['ccd_id_str'])

psf_params['ccd_cuts'] = ccd['ccd_cuts']
psf_params.remove_column('ccd_id_str')

psf_params.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr9/survey-ccds-decam-dr9-psfex.fits', overwrite=True)

print('All done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))

