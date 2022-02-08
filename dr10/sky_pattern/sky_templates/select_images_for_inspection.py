from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits
import healpy as hp
from astropy import wcs

from scipy.ndimage.filters import gaussian_filter
from pathlib import Path
from multiprocessing import Pool
import argparse

################################################################################

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'} 
plt.rcParams.update(params)


plot_dir = '/global/project/projectdirs/cosmo/www/temp/rongpu/dr10dev/compare_sky_corr'

skyrun = Table(fitsio.read('/global/cfs/cdirs/desi/users/schlafly/decals/skyrunsdr10-v3.fits'))
print('skyrun', len(skyrun))

mask = skyrun['filter']=='z'
skyrun = skyrun[mask]
print('skyrun', len(skyrun))

# Only plot exposures that are OK
mask = skyrun['ok']==True
skyrun = skyrun[mask]
print(len(skyrun))

for run in range(327, 430):

    print(run)

    mask = skyrun['run']==run
    band = skyrun['filter'][mask][0]

    np.random.seed(123+run)
    expnum_list = np.random.choice(skyrun['expnum'][mask], size=32, replace=False)

    expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)
    counter = 0
    for expnum in expnum_list:
        plot_path = os.path.join(plot_dir, band, '{}_{}_image_{}_medianscale.png'.format(band, run, expnum))
        if os.path.isfile(plot_path):
            os.system('cp {} /global/project/projectdirs/cosmo/www/temp/rongpu/dr10dev/compare_sky_corr/z/less/'.format(plot_path))
            os.system('cp {} /global/project/projectdirs/cosmo/www/temp/rongpu/dr10dev/compare_sky_corr/z/less/'.format(plot_path.replace('_medianscale', '_cp_original')))
            counter+=1
        if counter==4:
            break

    print('Done!!!!!!!!!!!!!!!!!!!!!')
