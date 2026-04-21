# Compare the approximate depth maps with the depth values in the tractor catalogs
# Also get the nominal efftime in my native ("approximate") units

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

# https://github.com/rongpu/Python/blob/master/user_modules/match_coord.py
sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']

columns = ['ra', 'dec', 'ebv', 'brick_primary', 'psfdepth_M411', 'psfdepth_M438', 'psfdepth_M464', 'psfdepth_M490', 'psfdepth_M517']

print('######################################## XMM-LSS ########################################')

cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_xmm.fits', columns=columns))
print(len(cat))

for col in cat.colnames:
    cat.rename_column(col, col.upper())

assert np.all(cat['BRICK_PRIMARY'])

ramin, ramax, decmin, decmax = cat['RA'].min()-0.05, cat['RA'].max()+0.05, cat['DEC'].min()-0.05, cat['DEC'].max()+0.05

for band in ibis_filters:
    cat[band+'_psfdepth'] = -2.5*(np.log10((5/np.sqrt(cat['PSFDEPTH_'+band])))-9)

median_zp_xmm = {}
for band in ibis_filters:
    tt = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/approximate_depth_maps/xmm-lss_{}.fits'.format(band)))
    mask = np.isfinite(tt['depth'])
    tt = tt[mask]

    idx1, idx2, d2d, d_ra, d_dec = match_coord(cat['RA'], cat['DEC'], tt['ra'], tt['dec'], search_radius=10., plot_q=False, verbose=False)
    # idx1, idx2, d2d, d_ra, d_dec = match_coord(ra1, dec1, ra2, dec2, search_radius=1., plot_q=True)

    x, y = cat[band+'_psfdepth'][idx1].copy(), tt['depth'][idx2].copy()
    median_zp_xmm[band] = np.median(y-x)
    print(band, median_zp_xmm[band])
print()

print('######################################## COSMOS ########################################')

cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_cosmos.fits', columns=columns))
print(len(cat))

for col in cat.colnames:
    cat.rename_column(col, col.upper())

assert np.all(cat['BRICK_PRIMARY'])

ramin, ramax, decmin, decmax = cat['RA'].min()-0.05, cat['RA'].max()+0.05, cat['DEC'].min()-0.05, cat['DEC'].max()+0.05

for band in ibis_filters:
    cat[band+'_psfdepth'] = -2.5*(np.log10((5/np.sqrt(cat['PSFDEPTH_'+band])))-9)

median_zp_cosmos = {}
for band in ibis_filters:
    tt = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/approximate_depth_maps/cosmos_{}.fits'.format(band)))
    mask = np.isfinite(tt['depth'])
    tt = tt[mask]

    idx1, idx2, d2d, d_ra, d_dec = match_coord(cat['RA'], cat['DEC'], tt['ra'], tt['dec'], search_radius=10., plot_q=False, verbose=False)
    # idx1, idx2, d2d, d_ra, d_dec = match_coord(ra1, dec1, ra2, dec2, search_radius=1., plot_q=True)

    x, y = cat[band+'_psfdepth'][idx1].copy(), tt['depth'][idx2].copy()
    median_zp_cosmos[band] = np.median(y-x)
    print(band, median_zp_cosmos[band])
print()

print('median_zp_xmm')
print(median_zp_xmm)
print()
print('median_zp_cosmos')
print(median_zp_cosmos)

median_zp = {}
print('XMM-LSS vs COSMOS:')
for band in ibis_filters:
    print('{}: {:.3f}'.format(band, median_zp_xmm[band] - median_zp_cosmos[band]))
    median_zp[band] = (median_zp_xmm[band]+median_zp_cosmos[band])/2
print()
print('Average ZP:')
print({k: float(f"{v:.3f}") for k, v in median_zp.items()})
print()

nominal_depth_with_zp = 25.0
nominal_efftime_dict = {}
for band in ibis_filters:
    nominal_efftime_dict[band] = 10**((nominal_depth_with_zp + median_zp[band]) / 2.5 * 2)
print('Nominal efftime:')
print({k: float(f"{v:.1f}") for k, v in nominal_efftime_dict.items()})
