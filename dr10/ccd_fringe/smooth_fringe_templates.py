from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool

from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

sys.path.append(os.path.expanduser('~/git/Python/useful/'))
from decam_postage_stamps import create_image


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

plot_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/plots'
fringe_templates_raw_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/fringe_templates_raw'
fringe_templates_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/fringe_templates'

vrange = 1e-3


def smooth_template(ccdnum):

    ccdname = ccdnamenumdict_inv[ccdnum]
    fn = os.path.join(fringe_templates_raw_dir, 'DECam_z_frg_{}_CCD{}_raw.fits'.format(ccdname, str(ccdnum).zfill(2)))
    print(fn)

    output_path = os.path.join(fringe_templates_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))

    if not os.path.isfile(fn):
        print(fn, 'does not exist')
        return None

    if os.path.isfile(output_path):
        print(output_path, 'already exists')
        return None

    img = fitsio.read(fn)

    # Determine how many NaN pixels on each edge
    a1, a2, b1, b2 = 14, 14, 14, 14
    while np.all(np.isnan(img[:a1+1])):
        a1 += 1
    while np.all(np.isnan(img[-a2-1:])):
        a2 += 1
    while np.all(np.isnan(img[:, :b1+1])):
        b1 += 1
    while np.all(np.isnan(img[:, -b2-1:])):
        b2 += 1
    if not (np.all(np.isnan(img[:a1])) and np.all(np.isnan(img[-a2:])) and np.all(np.isnan(img[:, :b1])) and np.all(np.isnan(img[:, -b2:]))):
        raise ValueError
    print('Edge pixels', a1, a2, b1, b2)

    # Trim the Nan Edges
    img = img[a1:-a2, b1:-b2]
    print(img.shape)

    sky_nmad = nmad(img)

    # 6-sigma masking
    mask = (img<-6*sky_nmad) | (img>6*sky_nmad)
    print('{} fraction of masked pixels: {:.5f}% ({})'.format(ccdname, 100*np.sum(mask)/np.sum(np.isfinite(mask)), np.sum(mask)))
    img[mask] = 0.

    # 3-sigma clipping
    mask = (img<-3*sky_nmad) | (img>3*sky_nmad)
    print('{} fraction of clipped pixels: {:.5f}% ({})'.format(ccdname, 100*np.sum(mask)/np.sum(np.isfinite(mask)), np.sum(mask)))
    mask = (img<-3*sky_nmad)
    img[mask] = -3*sky_nmad
    mask = (img>3*sky_nmad)
    img[mask] = 3*sky_nmad

    # Gaussian filtering
    kernel = Gaussian2DKernel(2)
    img_smooth = convolve(img, kernel, boundary='extend', nan_treatment='interpolate')
    img_smooth = np.pad(img_smooth, ((a1, a2), (b1, b2)), mode='constant', constant_values=0)

    mask = ~np.isfinite(img_smooth)
    img_smooth[mask] = 0

    fitsio.write(output_path, img_smooth, clobber=True)

    create_image((img_smooth).T, cmap='seismic', vmin=-vrange, vmax=vrange)
    plt.savefig(os.path.join(plot_dir, 'DECam_z_frg_{}_CCD{}.png'.format(ccdname, str(ccdnum).zfill(2))))
    plt.close()

    return None


n_processes = 64
with Pool(processes=n_processes) as pool:
    res = pool.map(smooth_template, np.arange(1, 63), chunksize=1)


