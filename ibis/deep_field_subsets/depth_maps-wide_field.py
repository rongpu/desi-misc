# Compute, save and plot approximate depth maps
# These depth maps are inspected to make sure that these approximate measurements are accurate enough for depth cuts in subsets

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import matplotlib.colors as colors
from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord, search_around

params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

field_names = ['COSMOS', 'XMM-LSS']
field_centers = {'COSMOS':[150.1, 2.182], 'XMM-LSS':[35.75, -4.75]}
field_size = 4.5  # degrees
n_points = int(2e5)

nominal_zp = {'M411': -21.117, 'M438': -21.308, 'M464': -21.381, 'M490': -21.431, 'M517': -21.446}
# nominal efftime for a depth of 25.0
nominal_efftime = {'M411': 1278.1, 'M438': 899.1, 'M464': 786.3, 'M490': 716.1, 'M517': 696.4}

ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']

half_width = 4094/2 * 0.262 / 3600
half_height = 2046/2 * 0.262 / 3600

radec_limits = [130, 225, -7.5, -2]
ramin, ramax, decmin, decmax = radec_limits

# area where the wide survey is complete
radec_limits_complete = [190, 215, -6.5, -4]
ramin1, ramax1, decmin1, decmax1 = radec_limits_complete

plot_dir = '/global/cfs/cdirs/cosmo/www/temp/rongpu/ibis/ibis_dr1_subsets/wide_field'

surveyccd_path = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/misc/survey-ccds-ibis-deep-and-wide-3-9.fits'
surveyccd_psfex_path = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/misc/survey-ccds-ibis-deep-and-wide-3-9-psfex.fits'

ccd = Table(fitsio.read(surveyccd_path))
ccd1 = Table(fitsio.read(surveyccd_psfex_path))
assert len(ccd)==len(ccd1) and np.all(ccd['expnum']==ccd1['expnum'])
for col in ccd1.colnames:
    if col in ccd.colnames:
        ccd1.remove_column(col)
ccd = hstack([ccd, ccd1])
print(len(ccd))

mask = (ccd['ccd_cuts']==0)
ccd = ccd[mask]
print(len(ccd))

ccd['median_ccdskycounts'], ccd['median_psf_fwhm'] = 0., 0.
tmp = Table()
tmp['expnum'], idx, tmp['count'] = np.unique(ccd['expnum'], return_counts=True, return_index=True)
for expnum in tmp['expnum']:
    mask = ccd['expnum']==expnum
    ccd['median_ccdskycounts'][mask] = np.median(ccd['ccdskycounts'][mask])
    ccd['median_psf_fwhm'][mask] = np.median(ccd['psf_fwhm'][mask])
ccd['efftime'] = 10**(0.4*ccd['zpt']-9) * ccd['exptime'] / (ccd['median_ccdskycounts'] * ccd['psf_fwhm']**2)

ccd_all = ccd.copy()


def get_ccd_mask(ccd_index):
    """Get mask of grid points covered by this CCD"""
    mask = (ra_grid > ccd['ra'][ccd_index] - half_width/np.cos(np.radians(dec_grid))) & \
           (ra_grid < ccd['ra'][ccd_index] + half_width/np.cos(np.radians(dec_grid))) & \
           (dec_grid > ccd['dec'][ccd_index] - half_height) & \
           (dec_grid < ccd['dec'][ccd_index] + half_height)
    return mask


def get_depth_map(band):

    # Make variables global for worker processes
    global ra_grid, dec_grid, ccd

    ccd = ccd_all.copy()
    mask = ccd['filter']==band
    ccd = ccd[mask]

    mask = (ccd['ra']>ramin) & (ccd['ra']<ramax) & (ccd['dec']>decmin) & (ccd['dec']<decmax)
    ccd = ccd[mask]

    def sample_ra_dec(ramin, ramax, decmin, decmax, n, seed=None):
        rng = np.random.default_rng(seed)
        ra = rng.uniform(ramin, ramax, n)
        # Uniform in solid angle requires uniform sampling in sin(dec)
        sin_min = np.sin(np.radians(decmin))
        sin_max = np.sin(np.radians(decmax))
        dec = np.degrees(np.arcsin(rng.uniform(sin_min, sin_max, n)))
        return ra, dec

    ra_grid, dec_grid = sample_ra_dec(ramin, ramax, decmin, decmax, n=n_points)

    efftime_grid = np.zeros(len(ra_grid))
    nexp_grid = np.zeros(len(ra_grid), dtype=int)
    seeing_ccd_grid = np.full((len(ccd), len(ra_grid)), np.nan)

    # Parallel mask computation
    with Pool(processes=128) as pool:
        ccd_mask_cube = np.array(pool.map(get_ccd_mask, range(len(ccd))))

    # Vectorized accumulation
    efftime_grid = np.dot(ccd_mask_cube.T, ccd['efftime'])
    nexp_grid = np.sum(ccd_mask_cube, axis=0)

    # Seeing grid (using mean instead of median for speed)
    for i in range(len(ccd)):
        seeing_ccd_grid[i, ccd_mask_cube[i]] = 0.262 * ccd['psf_fwhm'][i]
    seeing_grid = np.nanmean(seeing_ccd_grid, axis=0)
    depth_grid = 2.5*np.log10(np.sqrt(efftime_grid))

    mask_nan = nexp_grid==0
    efftime_grid[mask_nan] = np.nan
    seeing_grid[mask_nan] = np.nan
    depth_grid[mask_nan] = np.nan

    return nexp_grid, efftime_grid, depth_grid, seeing_grid, ra_grid, dec_grid


for band in ibis_filters:

    print(band)

    fig_efftime, axes_efftime = plt.subplots(1, 1, figsize=(20, 6))
    fig_efftimeh, axes_efftimeh = plt.subplots(1, 1, figsize=(7, 5))
    fig_nexp, axes_nexp = plt.subplots(1, 1, figsize=(20, 6))
    fig_nexph, axes_nexph = plt.subplots(1, 1, figsize=(7, 5))
    fig_depth, axes_depth = plt.subplots(1, 1, figsize=(20, 6))
    fig_depthh, axes_depthh = plt.subplots(1, 1, figsize=(7, 5))
    fig_seeing, axes_seeing = plt.subplots(1, 1, figsize=(20, 6))
    fig_seeingh, axes_seeingh = plt.subplots(1, 1, figsize=(7, 5))

    nexp_grid, efftime_grid, depth_grid, seeing_grid, ra, dec = get_depth_map(band)

    maps = Table()
    maps['ra'], maps['dec'] = ra, dec
    maps['nexp'], maps['efftime'], maps['depth'], maps['seeing'] = nexp_grid, efftime_grid, depth_grid, seeing_grid
    maps.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/approximate_depth_maps/wide_field_survey-ccds-3-9_{}.fits'.format(band), overwrite=True)

    # Get depths in actual zero points
    depth_grid -= nominal_zp[band]
    # Get efftimes in nominal efftimes
    efftime_grid /= nominal_efftime[band]

    mask_central = (ra>ramin1) & (ra<ramax1) & (dec>decmin1) & (dec<decmax1)

    ax = axes_efftime
    im1 = ax.scatter(ra, dec, c=efftime_grid, s=1.5,
                     cmap='jet', norm=colors.LogNorm())
    ax.add_patch(plt.Rectangle((ramin1, decmin1), ramax1-ramin1, decmax1-decmin1, edgecolor='red', facecolor='none'))
    ax.axis([ramax, ramin, decmin, decmax])
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')

    ax = axes_efftimeh
    mask = np.isfinite(efftime_grid)
    xx = efftime_grid
    ax.hist([xx[mask_central & mask], xx[mask]], 100, label=['complete area', 'all'], histtype='stepfilled', color=['C1', 'C0'])
    ax.set_xlabel('{}-band effective time'.format(band))
    ax.set_title('median={:.2f}; nmad={:.2f}'.format(np.median(xx[mask_central & mask]), nmad(xx[mask_central & mask])))
    ax.axvline(np.median(xx[mask_central & mask]), ls='--', lw=1, color='C3')
    ax.legend(loc='best')
    ax.grid(alpha=0.5)

    ax = axes_nexp
    im2 = ax.scatter(ra, dec, c=nexp_grid, s=1.5,
                     cmap='jet', norm=colors.LogNorm(vmin=1, vmax=200))
    ax.add_patch(plt.Rectangle((ramin1, decmin1), ramax1-ramin1, decmax1-decmin1, edgecolor='red', facecolor='none'))
    ax.axis([ramax, ramin, decmin, decmax])
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')

    ax = axes_nexph
    mask = np.isfinite(depth_grid)
    xx = nexp_grid
    bins = np.arange(xx.min()-0.5, xx.max()+1.5, 1)
    ax.hist([xx[mask_central & mask], xx[mask]], bins=bins, label=['complete area', 'all'], histtype='stepfilled', color=['C1', 'C0'])
    ax.set_xlabel('{}-band nexp'.format(band))
    ax.set_title('median={}; nmad={:.1f}'.format(np.median(xx[mask_central & mask]), nmad(xx[mask_central & mask])))
    ax.axvline(np.median(xx[mask_central & mask]), ls='--', lw=1, color='C3')
    ax.legend(loc='best')
    ax.grid(alpha=0.5)

    ax = axes_depth
    im4 = ax.scatter(ra, dec, c=depth_grid, s=1.5,
                     cmap='jet', vmin=23.5, vmax=26)
    ax.add_patch(plt.Rectangle((ramin1, decmin1), ramax1-ramin1, decmax1-decmin1, edgecolor='red', facecolor='none'))
    ax.axis([ramax, ramin, decmin, decmax])
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')

    ax = axes_depthh
    mask = np.isfinite(depth_grid)
    xx = depth_grid
    ax.hist([xx[mask_central & mask], xx[mask]], 100, label=['complete area', 'all'], histtype='stepfilled', color=['C1', 'C0'], range=(23.5, 26))
    ax.set_xlabel('{}-band depth (mag)'.format(band))
    ax.set_title('median={:.2f}; nmad={:.2f}'.format(np.median(xx[mask_central & mask]), nmad(xx[mask_central & mask])))
    ax.axvline(np.median(xx[mask_central & mask]), ls='--', lw=1, color='C3')
    ax.legend(loc='best')
    ax.grid(alpha=0.5)

    ax = axes_seeing
    im6 = ax.scatter(ra, dec, c=seeing_grid, s=1.5,
                     cmap='jet', vmin=0.8, vmax=1.6)
    ax.add_patch(plt.Rectangle((ramin1, decmin1), ramax1-ramin1, decmax1-decmin1, edgecolor='red', facecolor='none'))
    ax.axis([ramax, ramin, decmin, decmax])
    ax.set_xlabel('RA')
    ax.set_ylabel('DEC')

    ax = axes_seeingh
    mask = np.isfinite(seeing_grid)
    xx = seeing_grid
    ax.hist([xx[mask_central & mask], xx[mask]], 100, label=['complete area', 'all'], histtype='stepfilled', color=['C1', 'C0'])
    ax.set_xlabel('{}-band median seeing FWHM (arcsec)'.format(band))
    ax.set_title('median={:.2f}; nmad={:.2f}'.format(np.median(xx[mask_central & mask]), nmad(xx[mask_central & mask])))
    ax.axvline(np.median(xx[mask_central & mask]), ls='--', lw=1, color='C3')
    ax.legend(loc='best')
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
