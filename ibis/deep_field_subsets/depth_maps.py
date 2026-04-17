# Compute, save and plot crude depth maps
# These depth maps are inspected to make sure that these crude measurements are accurate enough for depth cuts in subsets

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import matplotlib.colors as colors

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord, search_around

params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

plot_dir = '/global/cfs/cdirs/cosmo/www/temp/rongpu/ibis/ibis_dr1_subsets/full_depth'

surveyccd_path = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/survey-ccds-ibis-dr1.fits'
surveyccd_psfex_path = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/survey-ccds-ibis-dr1-psfex.fits'

ccd = Table(fitsio.read(surveyccd_path))
ccd1 = Table(fitsio.read(surveyccd_psfex_path))
assert len(ccd)==len(ccd1) and np.all(ccd['expnum']==ccd1['expnum'])
for col in ccd1.colnames:
    if col in ccd.colnames:
        ccd1.remove_column(col)
ccd = hstack([ccd, ccd1])
print(len(ccd))

mask = (ccd['ccd_cuts']==0) & (ccd['failure'])
print(len(ccd))

ccd['efftime'] = 10**(0.4*ccd['zpt']-9) * ccd['exptime'] / (ccd['ccdskycounts'] * ccd['psf_fwhm']**2)

# ccd['median_ccdskycounts'], ccd['median_psf_fwhm'] = 0., 0.
# tmp = Table()
# tmp['expnum'], idx, tmp['count'] = np.unique(ccd['expnum'], return_counts=True, return_index=True)
# for expnum in tmp['expnum']:
#     mask = ccd['expnum']==expnum
#     ccd['median_ccdskycounts'][mask] = np.median(ccd['ccdskycounts'])
#     ccd['median_psf_fwhm'][mask] = np.median(ccd['psf_fwhm'])
# exp = ccd[idx]
# print(len(exp))
# exp = join(exp, tmp, keys='expnum')
# print(len(exp))

# exp['efftime'] = 10**(0.4*exp['zpt']-9) * exp['exptime'] / (exp['median_ccdskycounts'] * exp['median_psf_fwhm']**2)


field_names = ['COSMOS', 'XMM-LSS']
radec_limits = [[147.8, 152.5, -0.1, 4.5],
                [33.4, 38.1, -7, -2.4]]
radec_limits_new = []
for radec in radec_limits:
    ramin, ramax, decmin, decmax = radec
    radec_limits_new.append([ramin+0.2/np.cos(np.radians(decmin)), ramax-0.2/np.cos(np.radians(decmin)), decmin+0.2, decmax-0.2])
radec_limits = radec_limits_new

pointing_ra = np.array([150.1, 35.75])
pointing_dec = np.array([2.182, -4.75])

ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']

ccd_all = ccd.copy()

nominal_efftimes = {'M411': 200,
                        'M438': 200,
                        'M464': 200,
                        'M490': 200,
                        'M517': 200}

vrange_dict = {'M411': [100, 10000],
               'M438': [100, 10000],
               'M464': [100, 10000],
               'M490': [100, 10000],
               'M517': [100, 10000]}

half_width = 4094/2 * 0.262 / 3600
half_height = 2046/2 * 0.262 / 3600


def get_depth_map(band, field_index):

    ccd = ccd_all.copy()
    mask = ccd['filter']==band
    ccd = ccd[mask]
    radec = radec_limits[field_index]
    ramin, ramax, decmin, decmax = radec
    mask = (ccd['ra']>ramin) & (ccd['ra']<ramax) & (ccd['dec']>decmin) & (ccd['dec']<decmax)
    ccd = ccd[mask]

    ra_list = np.linspace(ramin, ramax, 400)
    dec_list = np.linspace(decmin, decmax, 400)
    # d_ra, d_dec = np.diff(ra_list)[0], np.diff(dec_list)[0]
    ra_grid, dec_grid = np.meshgrid(ra_list, dec_list)
    ra_grid = ra_grid.flatten()
    dec_grid = dec_grid.flatten()

    efftime_grid = np.zeros(len(ra_grid))
    nexp_grid = np.zeros(len(ra_grid))
    seeing_ccd_grid = np.full((len(ccd), len(ra_grid)), np.nan)

    for ccd_index in range(len(ccd)):
        mask = (ra_grid>ccd['ra'][ccd_index]-half_width/np.cos(np.radians(dec_grid))) & (ra_grid<ccd['ra'][ccd_index]+half_width/np.cos(np.radians(dec_grid))) \
             & (dec_grid>ccd['dec'][ccd_index]-half_height) & (dec_grid<ccd['dec'][ccd_index]+half_height)
        efftime_grid[mask] += ccd['efftime'][ccd_index]
        nexp_grid[mask]+=1
        seeing_ccd_grid[ccd_index, mask] = 0.262*ccd['psf_fwhm'][ccd_index]
    seeing_grid = np.nanmedian(seeing_ccd_grid, axis=0)

    depth_grid = 2.5*np.log10(np.sqrt(efftime_grid/nominal_efftimes[band]))

    mask_nan = nexp_grid==0
    nexp_grid[mask_nan] = np.nan
    efftime_grid[mask_nan] = np.nan
    seeing_grid[mask_nan] = np.nan
    depth_grid[mask_nan] = np.nan

    # Ignore pixels more than 2.05 deg from the field centers
    search_radius = 2.05*3600.
    idx1, idx2, d2d, _, _ = search_around(pointing_ra, pointing_dec, ra_grid, dec_grid, search_radius=search_radius, verbose=False)
    mask = np.full(len(efftime_grid), True)
    mask[idx2] = False
    efftime_grid[mask] = np.nan
    nexp_grid[mask] = np.nan
    seeing_grid[mask] = np.nan
    depth_grid[mask] = np.nan

    nexp_grid = nexp_grid.reshape(len(dec_list), len(ra_list))
    efftime_grid = efftime_grid.reshape(len(dec_list), len(ra_list))
    depth_grid = depth_grid.reshape(len(dec_list), len(ra_list))
    seeing_grid = seeing_grid.reshape(len(dec_list), len(ra_list))

    extent = np.array([ra_grid.min(), ra_grid.max(), dec_grid.min(), dec_grid.max()])

    return nexp_grid, efftime_grid, depth_grid, seeing_grid, extent, ra_grid, dec_grid


for band in ibis_filters:

    print(band)

    fig_efftime, axes_efftime = plt.subplots(1, 2, figsize=(12, 6))
    fig_efftimeh, axes_efftimeh = plt.subplots(1, 2, figsize=(9, 5))
    fig_nexp, axes_nexp = plt.subplots(1, 2, figsize=(12, 6))
    fig_nexph, axes_nexph = plt.subplots(1, 2, figsize=(9, 5))
    fig_depth, axes_depth = plt.subplots(1, 2, figsize=(12, 6))
    fig_depthh, axes_depthh = plt.subplots(1, 2, figsize=(9, 5))
    fig_seeing, axes_seeing = plt.subplots(1, 2, figsize=(12, 6))
    fig_seeingh, axes_seeingh = plt.subplots(1, 2, figsize=(9, 5))

    for field_index in range(len(field_names)):
        print('field', field_index)

        field_name = field_names[field_index]

        nexp_grid, efftime_grid, depth_grid, seeing_grid, extent, ra, dec = get_depth_map(band, field_index)

        maps = Table()
        maps['ra'], maps['dec'] = ra, dec
        maps['nexp'], maps['efftime'], maps['depth'], maps['seeing'] = nexp_grid.flatten(), efftime_grid.flatten(), depth_grid.flatten(), seeing_grid.flatten()
        maps.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/crude_depth_maps/{}_{}.fits'.format(field_name.lower(), band))

        ax = axes_efftime[field_index]
        im1 = ax.imshow(efftime_grid, extent=extent,
                   cmap='jet', norm=colors.LogNorm(vmin=vrange_dict[band][0], vmax=vrange_dict[band][1]), origin='lower', aspect='auto')
        ax.axis(extent[[1, 0, 2, 3]])
        ax.set_title(field_name)

        ax = axes_efftimeh[field_index]
        ax.hist(efftime_grid[np.isfinite(efftime_grid)].flatten(), 100)
        ax.set_xlabel('{}-band efftime'.format(band))
        ax.set_title(field_name)
        ax.grid(alpha=0.5)

        ax = axes_nexp[field_index]
        im2 = ax.imshow(nexp_grid, extent=extent,
                   cmap='jet', norm=colors.LogNorm(vmin=1, vmax=1e3), origin='lower', aspect='auto')
        ax.axis(extent[[1, 0, 2, 3]])
        ax.set_title(field_name)

        ax = axes_nexph[field_index]
        xx = nexp_grid[np.isfinite(depth_grid)].flatten()
        bins = np.arange(xx.min()-0.5, xx.max()+1.5, 1)
        ax.hist(xx, bins=bins)
        ax.set_xlabel('{}-band nexp'.format(band))
        ax.set_title(field_name)
        ax.grid(alpha=0.5)

        ax = axes_depth[field_index]
        im4 = ax.imshow(depth_grid, extent=extent,
                   cmap='jet', vmin=-0.2, vmax=2, origin='lower', aspect='auto')
        ax.axis(extent[[1, 0, 2, 3]])
        ax.set_title(field_name)

        ax = axes_depthh[field_index]
        ax.hist(depth_grid[np.isfinite(depth_grid)].flatten(), 100)
        ax.set_xlabel('{}-band depth with arbitrary ZP (mag)'.format(band))
        ax.set_title(field_name)
        ax.grid(alpha=0.5)

        ax = axes_seeing[field_index]
        im6 = ax.imshow(seeing_grid, extent=extent,
                   cmap='jet', vmin=0.8, vmax=1.6, origin='lower', aspect='auto')
        ax.axis(extent[[1, 0, 2, 3]])
        ax.set_title(field_name)

        ax = axes_seeingh[field_index]
        ax.hist(seeing_grid[np.isfinite(seeing_grid)].flatten(), 100)
        ax.set_xlabel('{}-band median seeing FWHM (arcsec)'.format(band))
        ax.set_title(field_name)
        ax.grid(alpha=0.5)

    fig_efftime.colorbar(im1, ax=axes_efftime, location='right', pad=0.015)
    fig_nexp.colorbar(im2, ax=axes_nexp, location='right', pad=0.015)
    fig_depth.colorbar(im4, ax=axes_depth, location='right', pad=0.015)
    fig_seeing.colorbar(im6, ax=axes_seeing, location='right', pad=0.015)

    fig_efftime.savefig(os.path.join(plot_dir, 'efftime_{}.png'.format(band)), bbox_inches='tight')
    fig_efftimeh.savefig(os.path.join(plot_dir, 'efftime_{}_hist.png'.format(band)), bbox_inches='tight')
    fig_nexp.savefig(os.path.join(plot_dir, 'nexp_{}.png'.format(band)), bbox_inches='tight')
    fig_nexph.savefig(os.path.join(plot_dir, 'nexp_{}_hist.png'.format(band)), bbox_inches='tight')
    fig_depth.savefig(os.path.join(plot_dir, 'depth_{}.png'.format(band)), bbox_inches='tight')
    fig_depthh.savefig(os.path.join(plot_dir, 'depth_{}_hist.png'.format(band)), bbox_inches='tight')
    fig_seeing.savefig(os.path.join(plot_dir, 'seeing_{}.png'.format(band)), bbox_inches='tight')
    fig_seeingh.savefig(os.path.join(plot_dir, 'seeing_{}_hist.png'.format(band)), bbox_inches='tight')
    plt.close()
