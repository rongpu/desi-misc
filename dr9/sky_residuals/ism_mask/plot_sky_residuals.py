from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp

sys.path.append(os.path.expanduser('~/git/desi-examples/imaging_systematics'))
from plot_healpix_map import plot_map

median_ranges = {'g': 0.003, 'r': 0.005, 'z': 0.007}
nmad_ranges = {'g': [0.001, 0.005], 'r': [0.001, 0.01], 'z': [0.003, 0.015]}

plot_dir = '/global/cfs/cdirs/desi/users/rongpu/imaging_sys/sky_residuals'

for nside in [256, 512]:

    npix = hp.nside2npix(nside)

    sky_south = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_{}_south.fits'.format(nside)))
    sky_south['PHOTSYS'] = 'S'
    sky_north = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_{}_north.fits'.format(nside)))
    sky_north['PHOTSYS'] = 'N'

    mask = (sky_north['DEC']>32.375)
    sky_north = sky_north[mask]
    mask = ~np.in1d(sky_south['HPXPIXEL'], sky_north['HPXPIXEL'])
    sky = vstack([sky_north, sky_south[mask]])

    nsource_min = 250 * (512/nside)**2
    mask = sky['nsource_g']>nsource_min
    mask &= sky['nsource_r']>nsource_min
    mask &= sky['nsource_z']>nsource_min
    print(np.sum(mask)/len(mask))
    sky = sky[mask]

    for band in ['g', 'r', 'z']:

        plot_path = os.path.join(plot_dir, 'dr9_sky_median_{}_{}.png'.format(band, nside))
        if not os.path.isfile(plot_path):
            plot_map(nside, sky['sky_median_'+band], pix=sky['HPXPIXEL'],
                     vmin=-median_ranges[band], vmax=median_ranges[band], cmap='seismic', nest=False,
                     title='sky_median_{} NSIDE={}'.format(band, nside), save_path=plot_path)

        plot_path = os.path.join(plot_dir, 'dr9_sky_nmad_{}_{}.png'.format(band, nside))
        if not os.path.isfile(plot_path):
            plot_map(nside, sky['sky_nmad_'+band], pix=sky['HPXPIXEL'],
                     vmin=nmad_ranges[band][0], vmax=nmad_ranges[band][1], cmap='viridis', nest=False,
                     title='sky_nmad_{} NSIDE={}'.format(band, nside), save_path=plot_path)
