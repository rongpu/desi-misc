from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


n_exposures = 100

deep_ra = np.array([150.1166])
deep_dec = np.array([2.2058])

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-defringed.fits'))
print(len(cat))

_, idx = np.unique(cat['expnum'], return_index=True)
exp = cat[idx].copy()
print(len(exp))

exp['exposure_efftime'] = 10**(0.4*exp['zpt']-9) * exp['exptime'] / (exp['median_ccdskycounts'] * exp['median_psf_fwhm']**2)

exp_search_radius = .25*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, exp['ra_bore'], exp['dec_bore'], search_radius=exp_search_radius, plot_q=False, keep_all_pairs=True)
exp = exp[idx2]
print(len(exp))

exp['keep'] = np.full(len(exp), False)

for band in ['g', 'r', 'i', 'z', 'Y']:
    mask = exp['filter']==band
    idx = np.where(mask)[0]
    idx1 = np.argsort(exp['exposure_efftime'][idx])
    exp['keep'][idx[idx1[-n_exposures:]]] = True

mask = exp['keep'].copy()
exp = exp[mask]
print(len(exp))

mask = np.in1d(cat['expnum'], exp['expnum'])
cat = cat[mask]
print('ccd', np.sum(mask))

cat.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/test/survey-ccds-dr10-deep-fields-v1-defringed-cosmos-{}.fits'.format(n_exposures), overwrite=True)
