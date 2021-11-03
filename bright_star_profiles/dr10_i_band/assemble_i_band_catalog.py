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


exp = Table(fitsio.read('/global/cscratch1/sd/rongpu/temp/dr10_i_band_exposures.fits'))

file_list = glob.glob('/global/cscratch1/sd/dstn/dr10pre/zpt/decam/DECam_CP-DR10c/*/*_ooi_i*-photom.fits')
print(len(file_list))

cat_stack = []
for fn_index, fn in enumerate(file_list):
    
    tmp = Table(fitsio.read(fn, columns=['phot_g_mean_mag']))
    mask = (tmp['phot_g_mean_mag']>16) & (tmp['phot_g_mean_mag']<18)
    idx = np.where(mask)[0]
    cat = Table(fitsio.read(fn, rows=idx))
    
    # annotated_fn = fn.replace('photom.fits', 'annotated.fits')
    # annotated = Table(fitsio.read())
    exp_index = np.where(exp['expnum']==cat['expnum'][0])[0][0]

    for col in ['airmass', 'fwhm', 'psfdepth', 'meansky', 'stdsky', 'ebv', 'psfdepth']:
        if col not in cat.colnames:
            cat[col] = exp[col][exp_index]

    cat_stack.append(cat)

    print(fn_index, len(cat))

cat = vstack(cat_stack)
cat.write('/global/cscratch1/sd/rongpu/temp/dr10_i_band_photom.fits')

###########################################################

cat = Table(fitsio.read('/global/cscratch1/sd/rongpu/temp/dr10_i_band_photom.fits'))
print(len(cat), len(np.unique(cat['gaia_sourceid'])))

# Quality cuts
mask = cat['airmass']<2
print(np.sum(~mask), np.sum(~mask)/len(mask))
mask &= cat['fwhm']<6
print(np.sum(~mask), np.sum(~mask)/len(mask))
mask &= cat['psfdepth']>22
print(np.sum(~mask), np.sum(~mask)/len(mask))
mask &= cat['meansky']<6
print(np.sum(~mask), np.sum(~mask)/len(mask))
mask &= cat['stdsky']<0.01
print(np.sum(~mask), np.sum(~mask)/len(mask))

mask &= cat['fracmasked']<0.1
print(np.sum(~mask), np.sum(~mask)/len(mask))

mask &= cat['psfmag']!=0
print(np.sum(~mask), np.sum(~mask)/len(mask))

mask &= np.isfinite(cat['phot_bp_mean_mag']) & np.isfinite(cat['phot_rp_mean_mag'])
print(np.sum(~mask), np.sum(~mask)/len(mask))

cat = cat[mask]
print(len(cat), len(np.unique(cat['gaia_sourceid'])))

# Remove duplidates keeping the higher EFFTIME objects

cat.sort('psfdepth', reverse=True)
_, idx_keep = np.unique(cat['gaia_sourceid'], return_index=True)

########### Make sure that the kept objects are the better ones ##########
mask_remove = ~np.in1d(np.arange(len(cat)), idx_keep)
cat_bad = cat[mask_remove].copy()
_, idx_bad = np.unique(cat_bad['gaia_sourceid'], return_index=True)
cat_bad = cat_bad[idx_bad]
cat_bad.sort('gaia_sourceid')

mask_better = np.in1d(cat['gaia_sourceid'][idx_keep], cat['gaia_sourceid'][mask_remove])
cat_better = cat[idx_keep][mask_better].copy()
cat_better.sort('gaia_sourceid')
print(np.all(cat_better['gaia_sourceid']==cat_bad['gaia_sourceid']))

plt.hist(cat_better['psfdepth'], 50, range=(22, 25), alpha=0.5)
plt.hist(cat_bad['psfdepth'], 50, range=(22, 25), alpha=0.5, color='r')
plt.show()

plt.hist(cat_better['psfdepth']-cat_bad['psfdepth'], 50, alpha=0.5)
plt.show()
print((cat_better['psfdepth']-cat_bad['psfdepth']).min())
############################################################################

cat = cat[idx_keep]
print(len(cat), len(np.unique(cat['gaia_sourceid'])))

columns = ['dpsfmag', 'psfmag', 'phot_g_mean_mag', 'phot_bp_mean_mag', 'phot_rp_mean_mag', 'ra_gaia', 'dec_gaia', 'gaia_sourceid', 'legacy_survey_mag', 'ra', 'dec', 'ebv']
cat = cat[columns]

# Add accurate per-object (not per-expopsure) EBV
import get_ebv_from_map
cat['ebv'] = get_ebv_from_map.get_ebv_from_map([cat['ra'], cat['dec']])

cat.write('/global/cscratch1/sd/rongpu/temp/dr10_i_band_photom-trim.fits', overwrite=True)

############################################################################

# Add parallel catalog with GAIA EDR3 photometry

cat = Table(fitsio.read('/global/cscratch1/sd/rongpu/temp/dr10_i_band_photom-trim.fits'))
print(len(cat))
# gaia = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_edr3_g_18_dr9.fits'))
gaia = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/useful/gaia_edr3_g_18.fits'))

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

cat.write('/global/cscratch1/sd/rongpu/temp/dr10_i_band_photom_gaia_edr3-trim.fits', overwrite=False)

