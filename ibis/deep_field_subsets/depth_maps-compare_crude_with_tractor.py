# Compare the approximate depth maps with the depth values in the tractor catalogs

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

###################################################### XMM-LSS ######################################################

cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_xmm.fits', columns=columns))
print(len(cat))

for col in cat.colnames:
    cat.rename_column(col, col.upper())

assert np.all(cat['BRICK_PRIMARY'])

ramin, ramax, decmin, decmax = cat['RA'].min()-0.05, cat['RA'].max()+0.05, cat['DEC'].min()-0.05, cat['DEC'].max()+0.05

for band in ibis_filters:
    cat[band+'_psfdepth'] = -2.5*(np.log10((5/np.sqrt(cat['PSFDEPTH_'+band])))-9)

mask = np.random.rand(len(cat))<0.2
for band in ibis_filters:
    plt.figure(figsize=(12, 10))
    plt.scatter(cat['RA'][mask], cat['DEC'][mask], c=cat[band+'_psfdepth'][mask], vmin=23.5, vmax=25.8, s=0.5, cmap='jet')
    plt.colorbar(aspect=30)
    plt.axis([ramax, ramin, decmin, decmax])
    plt.grid(alpha=0.5)
    plt.title(band+' 5-sigma PSF depth [mag]')
    plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/ibis/ibis_dr1_subsets/compare_approximate_maps_with_tractor/xmm_{}_tractor.png'.format(band))
    plt.close()

for band in ibis_filters:
    tt = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/approximate_depth_maps/xmm-lss_{}.fits'.format(band)))
    mask = np.isfinite(tt['depth'])
    tt = tt[mask]

    idx1, idx2, d2d, d_ra, d_dec = match_coord(cat['RA'], cat['DEC'], tt['ra'], tt['dec'], search_radius=10., plot_q=False)
    # idx1, idx2, d2d, d_ra, d_dec = match_coord(ra1, dec1, ra2, dec2, search_radius=1., plot_q=True)

    x, y = cat[band+'_psfdepth'][idx1].copy(), tt['depth'][idx2].copy()
    median_zp = np.median(y-x)
    print(median_zp)
    y -= median_zp
    # x -= np.median(x)
    # y -= np.median(y)
    plt.figure(figsize=(10, 10))
    plt.plot(x, y, '.', ms=0.5, alpha=0.1)
    plt.axis([24.5, 25.8, 24.5, 25.8])
    plt.grid(alpha=0.5)
    plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/ibis/ibis_dr1_subsets/compare_approximate_maps_with_tractor/xmm_{}_compare.png'.format(band))
    plt.close()

    mask = np.isfinite(tt['depth'])
    plt.figure(figsize=(12, 10))
    plt.scatter(tt['ra'][mask], tt['dec'][mask], c=tt['depth'][mask]-median_zp, vmin=23.5, vmax=25.8, s=0.5, cmap='jet')
    plt.colorbar(aspect=30)
    plt.axis([ramax, ramin, decmin, decmax])
    plt.grid(alpha=0.5)
    plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/ibis/ibis_dr1_subsets/compare_approximate_maps_with_tractor/xmm_{}_approximate.png'.format(band))
    plt.close()

###################################################### COSMOS ######################################################

cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_cosmos.fits', columns=columns))
print(len(cat))

for col in cat.colnames:
    cat.rename_column(col, col.upper())

assert np.all(cat['BRICK_PRIMARY'])

ramin, ramax, decmin, decmax = cat['RA'].min()-0.05, cat['RA'].max()+0.05, cat['DEC'].min()-0.05, cat['DEC'].max()+0.05

for band in ibis_filters:
    cat[band+'_psfdepth'] = -2.5*(np.log10((5/np.sqrt(cat['PSFDEPTH_'+band])))-9)

mask = np.random.rand(len(cat))<0.2
for band in ibis_filters:
    plt.figure(figsize=(12, 10))
    plt.scatter(cat['RA'][mask], cat['DEC'][mask], c=cat[band+'_psfdepth'][mask], vmin=23.5, vmax=25.8, s=0.5, cmap='jet')
    plt.colorbar(aspect=30)
    plt.axis([ramax, ramin, decmin, decmax])
    plt.grid(alpha=0.5)
    plt.title(band+' 5-sigma PSF depth [mag]')
    plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/ibis/ibis_dr1_subsets/compare_approximate_maps_with_tractor/cosmos_{}_tractor.png'.format(band))
    plt.close()

for band in ibis_filters:
    tt = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/approximate_depth_maps/cosmos_{}.fits'.format(band)))
    mask = np.isfinite(tt['depth'])
    tt = tt[mask]

    idx1, idx2, d2d, d_ra, d_dec = match_coord(cat['RA'], cat['DEC'], tt['ra'], tt['dec'], search_radius=10., plot_q=False)
    # idx1, idx2, d2d, d_ra, d_dec = match_coord(ra1, dec1, ra2, dec2, search_radius=1., plot_q=True)

    x, y = cat[band+'_psfdepth'][idx1].copy(), tt['depth'][idx2].copy()
    median_zp = np.median(y-x)
    print(median_zp)
    y -= median_zp
    # x -= np.median(x)
    # y -= np.median(y)
    plt.figure(figsize=(10, 10))
    plt.plot(x, y, '.', ms=0.5, alpha=0.1)
    plt.axis([24.5, 25.8, 24.5, 25.8])
    plt.grid(alpha=0.5)
    plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/ibis/ibis_dr1_subsets/compare_approximate_maps_with_tractor/cosmos_{}_compare.png'.format(band))
    plt.close()

    mask = np.isfinite(tt['depth'])
    plt.figure(figsize=(12, 10))
    plt.scatter(tt['ra'][mask], tt['dec'][mask], c=tt['depth'][mask]-median_zp, vmin=23.5, vmax=25.8, s=0.5, cmap='jet')
    plt.colorbar(aspect=30)
    plt.axis([ramax, ramin, decmin, decmax])
    plt.grid(alpha=0.5)
    plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/ibis/ibis_dr1_subsets/compare_approximate_maps_with_tractor/cosmos_{}_approximate.png'.format(band))
    plt.close()
