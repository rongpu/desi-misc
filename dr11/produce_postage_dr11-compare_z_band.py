from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool
sys.path.append(os.path.expanduser('~/git/Python/useful/'))
from decam_postage_stamps import decam_postage_stamp


def wrapper(index):
    img_path = os.path.join(img_dir, exp_dr11['image_filename'][index])
    band = exp_dr11['filter'][index]
    vrange = image_vrange[band] * (exp_dr11['exptime'][index] * exp_dr11['ccdskycounts'][index] / 3000.)
    plot_path = os.path.join(plot_dir, str(int(exp_dr11['mjd_obs'][index]))+'_'+os.path.basename(exp_dr11['image_filename'][index]).replace('.fits.fz', '_{}_{}_dr11.png'.format(exp_dr11['expnum'][index], exp_dr11['plver'][index])))
    decam_postage_stamp(img_path, binsize=60, plot_path=plot_path, save_path=None, vrange=[-vrange, vrange], dr8=False, median=True,
        blob_mask=False, ood_mask=True, show=False, cmap='gray')

def wrapper1(index):
    img_path = os.path.join(img_dir, exp_dr9['image_filename'][index])
    band = exp_dr9['filter'][index]
    vrange = image_vrange[band] * (exp_dr9['exptime'][index] * exp_dr9['ccdskycounts'][index] / 3000.)
    plot_path = os.path.join(plot_dir, str(int(exp_dr9['mjd_obs'][index]))+'_'+os.path.basename(exp_dr9['image_filename'][index]).replace('.fits.fz', '_{}_{}_dr9.png'.format(exp_dr9['expnum'][index], exp_dr9['plver'][index])))
    decam_postage_stamp(img_path, binsize=60, plot_path=plot_path, save_path=None, vrange=[-vrange, vrange], dr8=False, median=True,
        blob_mask=False, ood_mask=True, show=False, cmap='gray')


n_processes = 128

image_vrange = {'u':5/2, 'g':5/2, 'r':6/2, 'i':10/2, 'z':30/2, 'Y':30/2}

# band = 'i'
# band = str(sys.argv[1])

img_dir = '/dvs_ro/cfs/cdirs/cosmo/staging'

plot_dir = '/global/cfs/cdirs/cosmo/www/temp/rongpu/dr11/postage_stamps_compare_with_dr9/'

exp_dr11 = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/dr11_postage/dr11_z_for_inspection.fits'))
exp_dr9 = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/dr11_postage/dr9_z_for_inspection.fits'))

np.random.seed(99)
idx = np.random.choice(len(exp_dr11), size=32, replace=False)
exp_dr11 = exp_dr11[idx]
exp_dr9 = exp_dr9[idx]

with Pool(processes=n_processes) as pool:
    res = pool.map(wrapper, np.arange(len(exp_dr11)))

with Pool(processes=n_processes) as pool:
    res = pool.map(wrapper1, np.arange(len(exp_dr9)))
