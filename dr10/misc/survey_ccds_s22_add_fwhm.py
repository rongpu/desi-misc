# Add PSFEx FWHM and Moffat fits

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool


psfex_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/calib/patched-psfex'

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-s22.fits'))


def get_psfex_params(index):
    psfex_fn = cat['image_filename'][index].replace('.fits.fz', '-psfex.fits')
    psfex_fn = os.path.join(psfex_dir, psfex_fn)
    if not os.path.isfile(psfex_fn):
        return None
    tmp = Table(fitsio.read(psfex_fn))
    if not 'moffat_alpha' in tmp.colnames:
        return None
    mask = tmp['ccdname']=='S22'
    if np.sum(mask)==0:
        return None
    tmp = tmp[mask]
    tmp = tmp[['expnum', 'psf_fwhm', 'moffat_alpha', 'moffat_beta', 'failure']]
    return tmp


print('Start!')

n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(get_psfex_params, np.arange(len(cat)))

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

psf_params = vstack(res)
psf_params.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-s22-psfex-moffat-params.fits')
