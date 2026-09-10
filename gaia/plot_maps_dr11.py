# source /global/common/software/desi/desi_environment.sh 22.2

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/desi-examples/imaging_systematics'))
from plot_healpix_map import plot_map

params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

## DR11 South

dr11_south_offsets = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_xp_dr11_south_offset_maps_256.fits'))

for band in ['g', 'r', 'i', 'z']:
    plot_map(256, np.array(dr11_south_offsets[band+'mag_n_objects']).astype(float), pix=dr11_south_offsets['HPXPIXEL'], dpi=400, xsize=6000, cmap='gray_r',
         save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_south_density_{}_256.png'.format(band), vmin=0, vmax=200,
         cbar_label='# stars in pixel ({} band)'.format(band))

for band in ['g', 'r', 'i', 'z']:
    mask = np.isfinite(dr11_south_offsets[band+'mag_diff_median'])
    plot_map(256, 1000*dr11_south_offsets[band+'mag_diff_median'][mask], pix=dr11_south_offsets['HPXPIXEL'][mask], dpi=400, xsize=6000, cmap='seismic',
             vmin=-30, vmax=30, cbar_label='${}_\mathrm{{LS}}-{}_\mathrm{{Gaia,synth}}$ (mmag)'.format(band, band),
             save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_south_median_offset_{}_256.png'.format(band))

for band in ['g', 'r', 'i', 'z']:
    mask = np.isfinite(dr11_south_offsets[band+'mag_diff_median'])
    plot_map(256, 1000*dr11_south_offsets[band+'mag_diff_median'][mask], pix=dr11_south_offsets['HPXPIXEL'][mask], dpi=400, xsize=6000, cmap='gray',
             vmin=-30, vmax=30, cbar_label='${}_\mathrm{{LS}}-{}_\mathrm{{Gaia,synth}}$ (mmag)'.format(band, band),
             save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_south_median_offset_{}_256_gray.png'.format(band))

## DR11 North

dr11_north_offsets = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_xp_dr11_north_offset_maps_256.fits'))

for band in ['g', 'r', 'z']:
    plot_map(256, np.array(dr11_north_offsets[band+'mag_n_objects']).astype(float), pix=dr11_north_offsets['HPXPIXEL'], dpi=400, xsize=6000, cmap='gray_r',
         save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_north_density_{}_256.png'.format(band), vmin=0, vmax=200,
         cbar_label='# stars in pixel ({} band)'.format(band))

for band in ['g', 'r', 'z']:
    mask = np.isfinite(dr11_north_offsets[band+'mag_diff_median'])
    plot_map(256, 1000*dr11_north_offsets[band+'mag_diff_median'][mask], pix=dr11_north_offsets['HPXPIXEL'][mask], dpi=400, xsize=6000, cmap='seismic',
             vmin=-30, vmax=30, cbar_label='${}_\mathrm{{LS}}-{}_\mathrm{{Gaia,synth}}$ (mmag)'.format(band, band),
             save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_north_median_offset_{}_256.png'.format(band))

for band in ['g', 'r', 'z']:
    mask = np.isfinite(dr11_north_offsets[band+'mag_diff_median'])
    plot_map(256, 1000*dr11_north_offsets[band+'mag_diff_median'][mask], pix=dr11_north_offsets['HPXPIXEL'][mask], dpi=400, xsize=6000, cmap='gray',
             vmin=-30, vmax=30, cbar_label='${}_\mathrm{{LS}}-{}_\mathrm{{Gaia,synth}}$ (mmag)'.format(band, band),
             save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_north_median_offset_{}_256_gray.png'.format(band))

## DR11 North+South

dr11_south_offsets = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_xp_dr11_south_offset_maps_256.fits'))
dr11_north_offsets = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_xp_dr11_north_offset_maps_256.fits'))

mask = (dr11_north_offsets['DEC']>32.375) & (dr11_north_offsets['RA']>90) & (dr11_north_offsets['RA']<300)
dr11_north_offsets = dr11_north_offsets[mask]

mask = ~np.in1d(dr11_south_offsets['HPXPIXEL'], dr11_north_offsets['HPXPIXEL'])
dr11_south_offsets = dr11_south_offsets[mask]

dr11_offsets = vstack([dr11_south_offsets, dr11_north_offsets])
dr11_offsets.sort('HPXPIXEL')

for band in ['g', 'r', 'z']:
    plot_map(256, np.array(dr11_offsets[band+'mag_n_objects']).astype(float), pix=dr11_offsets['HPXPIXEL'], dpi=400, xsize=6000, cmap='gray_r',
         save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_density_{}_256.png'.format(band), vmin=0, vmax=200,
         cbar_label='# stars in pixel ({} band)'.format(band))

for band in ['g', 'r', 'z']:
    mask = np.isfinite(dr11_offsets[band+'mag_diff_median'])
    plot_map(256, 1000*dr11_offsets[band+'mag_diff_median'][mask], pix=dr11_offsets['HPXPIXEL'][mask], dpi=400, xsize=6000, cmap='seismic',
             vmin=-30, vmax=30, cbar_label='${}_\mathrm{{LS}}-{}_\mathrm{{Gaia,synth}}$ (mmag)'.format(band, band),
             save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_median_offset_{}_256.png'.format(band))

for band in ['g', 'r', 'z']:
    mask = np.isfinite(dr11_offsets[band+'mag_diff_median'])
    plot_map(256, 1000*dr11_offsets[band+'mag_diff_median'][mask], pix=dr11_offsets['HPXPIXEL'][mask], dpi=400, xsize=6000, cmap='gray',
             vmin=-30, vmax=30, cbar_label='${}_\mathrm{{LS}}-{}_\mathrm{{Gaia,synth}}$ (mmag)'.format(band, band),
             save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_median_offset_{}_256_gray.png'.format(band))

for band in ['g', 'r', 'z']:
    mask = np.isfinite(dr11_offsets[band+'mag_diff_mean'])
    plot_map(256, 1000*dr11_offsets[band+'mag_diff_mean'][mask], pix=dr11_offsets['HPXPIXEL'][mask], dpi=400, xsize=6000, cmap='seismic',
             vmin=-30, vmax=30, cbar_label='${}_\mathrm{{LS}}-{}_\mathrm{{Gaia,synth}}$ (mmag)'.format(band, band),
             save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_mean_offset_{}_256.png'.format(band))

for band in ['g', 'r', 'z']:
    mask = np.isfinite(dr11_offsets[band+'mag_diff_mean'])
    plot_map(256, 1000*dr11_offsets[band+'mag_diff_mean'][mask], pix=dr11_offsets['HPXPIXEL'][mask], dpi=400, xsize=6000, cmap='gray',
             vmin=-30, vmax=30, cbar_label='${}_\mathrm{{LS}}-{}_\mathrm{{Gaia,synth}}$ (mmag)'.format(band, band),
             save_path='/global/cfs/cdirs/cosmo/www/temp/rongpu/gaia/ls_comparison/dr11/gaia_xp_dr11_mean_offset_{}_256_gray.png'.format(band))

