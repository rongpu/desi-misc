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

dr10_offsets = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_xp_dr10_offset_maps_256.fits'))

for band in ['g', 'r', 'i', 'z']:
    plot_map(256, dr10_offsets['HPXPIXEL'], np.array(dr10_offsets[band+'mag_n_objects']).astype(float), dpi=1200, xsize=6000, cmap='gray_r',
             save_path='plots/gaia_xp_dr10_density_{}_256.png'.format(band), vmin=0, vmax=200)

for band in ['g', 'r', 'i', 'z']:
    mask = np.isfinite(dr10_offsets[band+'mag_diff'])
    plot_map(256, dr10_offsets['HPXPIXEL'][mask], dr10_offsets[band+'mag_diff'][mask], dpi=1200, xsize=6000, cmap='seismic',
             vmin=-0.04, vmax=0.04,
             save_path='plots/gaia_xp_dr10_offset_{}_256.png'.format(band))

for band in ['g', 'r', 'i', 'z']:
    mask = np.isfinite(dr10_offsets[band+'mag_diff'])
    plot_map(256, dr10_offsets['HPXPIXEL'][mask], dr10_offsets[band+'mag_diff'][mask], dpi=1200, xsize=6000, cmap='gray',
             vmin=-0.03, vmax=0.03,
             save_path='plots/gaia_xp_dr10_offset_{}_256_gray.png'.format(band))
