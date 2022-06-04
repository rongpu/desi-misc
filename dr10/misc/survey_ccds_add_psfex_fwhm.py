# Add PSFEx FWHM values

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool


psfex_unpatched_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/calib/unpatched-psfex'
psfex_patched_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/calib/patched-psfex'

# Unique exposures
cat = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', columns=['expnum', 'image_filename', 'ccdname']))
_, idx = np.unique(cat['expnum'], return_index=True)
exp = cat[idx].copy()

cat['ccd_id_str'] = np.char.add(np.array(cat['expnum']).astype(str), cat['ccdname'])

def get_psfex_params(index):
    psfex_filename = exp['image_filename'][index].replace('.fits.fz', '-psfex.fits')
    psfex_fn = os.path.join(psfex_patched_dir, psfex_filename)
    if not os.path.isfile(psfex_fn):
        psfex_fn = os.path.join(psfex_unpatched_dir, psfex_filename)
        if not os.path.isfile(psfex_fn):
            return None
    tt = Table(fitsio.read(psfex_fn))
    if len(tt)==0:
        return None
    tt['ccd_id_str'] = np.char.add(np.array(tt['expnum']).astype(str), tt['ccdname'])
    tt['median_psf_fwhm'] = np.median(tt['psf_fwhm'])
    if 'moffat_alpha' in tt.colnames:
        tt = tt[['ccd_id_str', 'psf_fwhm', 'median_psf_fwhm', 'moffat_alpha', 'moffat_beta', 'failure']]
    else:
        tt = tt[['ccd_id_str', 'psf_fwhm', 'median_psf_fwhm']]
    return tt

print('Start!')

n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(get_psfex_params, np.arange(len(exp)))

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

psf_params = vstack(res).filled(0)
cat1 = cat[['expnum', 'ccdname', 'ccd_id_str']]
psf_params = join(cat1, psf_params, keys='ccd_id_str', join_type='left')
psf_params.remove_column('ccd_id_str')

psf_params.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-psfex-fwhm.fits', overwrite=True)
