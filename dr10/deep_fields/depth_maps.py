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

surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1.fits'
plot_dir = '/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/new/all_ccds/'

ccd = Table(fitsio.read(surveyccd_path))
print(len(ccd))

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx]

exp['efftime'] = 10**(0.4*exp['zpt']-9) * exp['exptime'] / (exp['median_ccdskycounts'] * exp['median_psf_fwhm']**2)

field_names = ['COSMOS', 'S1 & S2', 'X1 & X2 & X3 (XMM-LSS)', 'C1 & C2 & C3', 'E1 & E2']
radec_limits = [[147.8, 152.5, -0.1, 4.5],
                [38.9, 45.1, -3.3, 2.3],
                [32.2, 38.8, -8.7, -2.3],
                [50, 56.9, -31.5, -24.8],
                [4.7, 12.7, -46.3, -40.7]]
radec_limits_new = []
for radec in radec_limits:
    ramin, ramax, decmin, decmax = radec
    radec_limits_new.append([ramin+0.2/np.cos(np.radians(decmin)), ramax-0.2/np.cos(np.radians(decmin)), decmin+0.2, decmax-0.2])
radec_limits = radec_limits_new

pointing_ra = np.array([150.1166, 54.2743, 54.2743, 52.6484, 34.4757, 35.6645, 36.4500, 42.8200, 41.1944, 7.8744, 9.5000])
pointing_dec = np.array([2.2058, -27.1116, -29.0884, -28.1000, -4.9295, -6.4121, -4.6000, 0.0000, -0.9884, -43.0096, -43.9980])
pointing_names = np.array(['COSMOS', 'C1', 'C2', 'C3', 'X1', 'X2', 'X3', 'S1', 'S2', 'E1', 'E2'])

ccd_all = ccd.copy()

des_nominal_efftimes = {'g': 239.6,
                        'r': 139.9,
                        'i': 53.4,
                        'z': 19.0,
                        'Y': 4.2}

vrange_dict = {'g': [100, 10000],
               'r': [100, 4000],
               'i': [50, 2000],
               'z': [20, 2000],
               'Y': [1, 500]}

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

    depth_grid = 2.5*np.log10(np.sqrt(efftime_grid/des_nominal_efftimes[band]))

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

    return nexp_grid, efftime_grid, depth_grid, seeing_grid, extent


for band in ['g', 'r', 'i', 'z', 'Y']:

    print(band)

    fig_efftime, axes_efftime = plt.subplots(1, 5, figsize=(42, 6))
    fig_efftimeh, axes_efftimeh = plt.subplots(1, 5, figsize=(30, 5))
    fig_nexp, axes_nexp = plt.subplots(1, 5, figsize=(42, 6))
    fig_nexph, axes_nexph = plt.subplots(1, 5, figsize=(30, 5))
    fig_depth, axes_depth = plt.subplots(1, 5, figsize=(42, 6))
    fig_depthh, axes_depthh = plt.subplots(1, 5, figsize=(30, 5))
    fig_seeing, axes_seeing = plt.subplots(1, 5, figsize=(42, 6))
    fig_seeingh, axes_seeingh = plt.subplots(1, 5, figsize=(30, 5))

    for field_index in range(len(field_names)):
        print('field', field_index)

        field_name = field_names[field_index]

        nexp_grid, efftime_grid, depth_grid, seeing_grid, extent = get_depth_map(band, field_index)

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
        ax.hist(nexp_grid[np.isfinite(depth_grid)].flatten(), 100)
        ax.set_xlabel('{}-band nexp'.format(band))
        ax.set_title(field_name)
        ax.grid(alpha=0.5)

        ax = axes_depth[field_index]
        im4 = ax.imshow(depth_grid, extent=extent,
                   cmap='jet', vmin=0, vmax=3, origin='lower', aspect='auto')
        ax.axis(extent[[1, 0, 2, 3]])
        ax.set_title(field_name)

        ax = axes_depthh[field_index]
        ax.hist(depth_grid[np.isfinite(depth_grid)].flatten(), 100)
        ax.set_xlabel('{}-band depth relative to DES (mag)'.format(band))
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
