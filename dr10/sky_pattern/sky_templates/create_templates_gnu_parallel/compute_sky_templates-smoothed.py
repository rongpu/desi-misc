from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter

from multiprocessing import Pool
import argparse
from pathlib import Path

from scipy.ndimage.filters import gaussian_filter

from scipy.interpolate import interp2d

params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

n_processes = 32
diagnostic_touch = True


nmad = lambda x: 1.4826 * np.nanmedian(np.abs(x-np.nanmedian(x)))

ccdnamenumdict = {'S1': 25, 'S2': 26, 'S3': 27, 'S4': 28,
                  'S5': 29, 'S6': 30, 'S7': 31,
                  'S8': 19, 'S9': 20, 'S10': 21, 'S11': 22, 'S12': 23,
                  'S13': 24,
                  'S14': 13, 'S15': 14, 'S16': 15, 'S17': 16, 'S18': 17,
                  'S19': 18,
                  'S20': 8, 'S21': 9, 'S22': 10, 'S23': 11, 'S24': 12,
                  'S25': 4, 'S26': 5, 'S27': 6, 'S28': 7,
                  'S29': 1, 'S30': 2, 'S31': 3,
                  'N1': 32, 'N2': 33, 'N3': 34, 'N4': 35,
                  'N5': 36, 'N6': 37, 'N7': 38,
                  'N8': 39, 'N9': 40, 'N10': 41, 'N11': 42, 'N12': 43,
                  'N13': 44,
                  'N14': 45, 'N15': 46, 'N16': 47, 'N17': 48, 'N18': 49,
                  'N19': 50,
                  'N20': 51, 'N21': 52, 'N22': 53, 'N23': 54, 'N24': 55,
                  'N25': 56, 'N26': 57, 'N27': 58, 'N28': 59,
                  'N29': 60, 'N30': 61, 'N31': 62,
                  }
ccdnamenumdict_inv = {aa: bb for bb, aa in ccdnamenumdict.items()}

ccdnum_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
               18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
               35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
               52, 53, 54, 55, 56, 57, 58, 59, 60, 62]

binsize = 4  # 4x4 downsizing


smoothing_scale = 30

runs_without_n15 = [392, 393, 394]  # hot spot
runs_with_cp_fringe = [418, 419, 420, 421, 422]

templates_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_templates'

skyrun = Table.read('/global/cfs/cdirs/desi/users/schlafly/decals/skyrunsdr10-v3.fits')


def compute_smooth_sky(run):

    # if run in runs_with_cp_fringe:
    if run not in runs_with_cp_fringe:
        return None

    mask = skyrun['run']==run
    band = skyrun['filter'][mask][0]

    raw_path = os.path.join(templates_dir, 'sky_raw_{}_{}.fits'.format(band, run))
    output_path = os.path.join(templates_dir, 'sky_template_{}_{}.fits'.format(band, run))

    if not os.path.isfile(raw_path):
        return None

    if os.path.isfile(output_path):
        return None

    print('run', run)

    hdul_template = fitsio.FITS(output_path, mode='rw', clobber=True)
    hdul_template.write(data=None)  # first HDU is empty

    for ccdnum in ccdnum_list:

        ccdname = ccdnamenumdict_inv[ccdnum]

        if ccdname=='S7' or ccdname=='S30':
            continue

        if ccdname=='N15' and run in runs_without_n15:
            continue

        try:
            img = fitsio.read(raw_path, ext=ccdname)
        except:
            continue

        img = img[3:-3, 3:-3]   # remove NaN edges

        mask = ~np.isfinite(img)
        img[mask] = 0.

        # 5-sigma clipping
        sky_nmad = nmad(img)  # sky level
        mask = (img<-5*sky_nmad) | (img>5*sky_nmad)
        print(ccdname, '5-sigma clipped pixels: {} ({:.2f}%)'.format(np.sum(mask), np.sum(mask)/len(mask.flatten())*100))
        # img[mask] = 0
        mask = (img<-5*sky_nmad)
        img[mask] = -5*sky_nmad
        mask = (img>5*sky_nmad)
        img[mask] = 5*sky_nmad

        naxis1, naxis2 = 2046, 4094
        pix_x_orig, pix_y_orig = np.arange(naxis1), np.arange(naxis2)
        pix_x, pix_y = np.arange(-1, naxis1+1), np.arange(-1, naxis2+1)
        pix_x = np.mean(pix_x.reshape(len(pix_x)//binsize, -1), axis=1)
        pix_y = np.mean(pix_y.reshape(len(pix_y)//binsize, -1), axis=1)
        pix_x = pix_x[3:-3]
        pix_y = pix_y[3:-3]

        img_smooth = gaussian_filter(img, smoothing_scale, mode='reflect')

        f = interp2d(pix_y, pix_x, img_smooth.flatten(), kind='linear')
        img_smooth_orig = f(pix_y_orig, pix_x_orig).T

        hdul_template.write(data=img_smooth_orig, extname=ccdname, compress='rice')

    hdul_template.close()


run_list = np.arange(429+1)
with Pool(processes=n_processes) as pool:
    res = pool.map(compute_smooth_sky, run_list, chunksize=1)

