from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp

sys.path.append(os.path.expanduser('~/git/Python/'))
import get_ebv_from_map

############################################################################

ccd = Table(fitsio.read('/global/cscratch1/sd/dstn/survey-ccds-Y.fits'))

# Quality cuts
mask = ccd['ccdzpt']!=0
print(np.sum(mask)/len(mask))
ccd = ccd[mask]

mask = ccd['airmass']<1.5
print(np.sum(mask), np.sum(mask)/len(mask))
mask &= ccd['fwhm']<6
print(np.sum(mask), np.sum(mask)/len(mask))
ccd = ccd[mask]

mask = ccd['ccdskycounts']<30
print(np.sum(mask)/len(mask))
ccd = ccd[mask]

mask = ccd['skyrms']<0.4
print(np.sum(mask)/len(mask))
ccd = ccd[mask]

mask = ccd['ccdzpt']>23.5
print(np.sum(mask)/len(mask))
ccd = ccd[mask]

# keep one entry per exposure
mask = ccd['ccdname']=='N10'
print(np.sum(mask)/len(mask))
ccd = ccd[mask]
print(len(ccd))

# Down-select 500 exposures
np.random.seed(63821)
idx = np.random.choice(len(ccd), size=500, replace=False)
ccd = ccd[idx]
print(len(ccd))

############################################################################

cat_stack = []
for index in range(len(ccd)):
    fn = ccd['image_filename'][index].replace('.fits.fz', '-photom.fits')
    fn = os.path.join('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/zpt/', fn)
    tmp = Table(fitsio.read(fn, columns=['phot_g_mean_mag']))
    mask = (tmp['phot_g_mean_mag']>16) & (tmp['phot_g_mean_mag']<18)
    idx = np.where(mask)[0]
    cat = Table(fitsio.read(fn, rows=idx))
    cat_stack.append(cat)

cat = vstack(cat_stack)

mask = cat['fracmasked']<0.1
print(np.sum(mask)/len(mask))
cat = cat[mask]

mask = cat['psfmag']!=0
print(np.sum(mask)/len(mask))
cat = cat[mask]

print(len(cat), len(np.unique(cat['gaia_sourceid'])))

# Remove duplidates keeping objects with smaller magnitude errors
cat.sort('dpsfmag', reverse=False)
_, idx_keep = np.unique(cat['gaia_sourceid'], return_index=True)

cat = cat[idx_keep]
print(len(cat), len(np.unique(cat['gaia_sourceid'])))

# Add accurate per-object (not per-expopsure) EBV
cat['ebv'] = get_ebv_from_map.get_ebv_from_map([cat['ra'], cat['dec']])

columns = ['dpsfmag', 'psfmag', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'ra_gaia', 'dec_gaia', 'gaia_sourceid', 'legacy_survey_mag', 'ra', 'dec', 'ebv']
cat = cat[columns]

cat.write('/global/cscratch1/sd/rongpu/temp/dr10_Y_band_photom_gaia_dr2-trim.fits', overwrite=False)

############################################################################

# Add parallel catalog with GAIA EDR3 photometry

cat = Table(fitsio.read('/global/cscratch1/sd/rongpu/temp/dr10_Y_band_photom_gaia_dr2-trim.fits'))
gaia = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_edr3_g_18_dr9.fits'))

nside = 128
cat_idx = hp.ang2pix(nside, cat['ra_gaia'], cat['dec_gaia'], lonlat=True)
gaia_idx = hp.ang2pix(nside, gaia['RA'], gaia['DEC'], lonlat=True)
mask = np.in1d(gaia_idx, cat_idx)

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

idx1, idx2, d2d, d_ra, d_dec = match_coord(cat['ra_gaia'], cat['dec_gaia'], gaia['RA'], gaia['DEC'], search_radius=0.1, plot_q=False)
cat = cat[idx1]
gaia = gaia[idx2]

cat.remove_columns(['phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'ra_gaia', 'dec_gaia', 'gaia_sourceid'])
gaia.rename_columns(['RA', 'DEC'], ['RA_GAIA', 'DEC_GAIA'])
cat = hstack([cat, gaia])

cat.write('/global/cscratch1/sd/rongpu/temp/dr10_Y_band_photom_gaia_edr3-trim.fits', overwrite=False)
