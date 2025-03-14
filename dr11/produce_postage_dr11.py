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


n_processes = 128

image_vrange = {'u':5, 'g':5, 'r':6, 'i':10, 'z':30, 'Y':30}

# band = 'i'
# band = str(sys.argv[1])

img_dir = '/dvs_ro/cfs/cdirs/cosmo/staging'

plot_dir = '/global/cfs/cdirs/cosmo/www/temp/rongpu/dr11/postage_stamps/'

exp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/tmp/dr11_exposures_for_inspection1.fits'))

np.random.seed(99)
idx = np.random.choice(len(exp), size=128, replace=False)
exp = exp[idx]

exp.sort('expnum')

def wrapper(index):
    img_path = os.path.join(img_dir, exp['image_filename'][index])
    band = exp['filter'][index]
    vrange = image_vrange[band] * (exp['exptime'][index] * exp['ccdskycounts'][index] / 3000.)
    plot_path = os.path.join(plot_dir, str(int(exp['mjd_obs'][index]))+'_'+exp['plver'][index]+'_'+os.path.basename(exp['image_filename'][index]).replace('.fits.fz', '_{}.png'.format(exp['expnum'][index])))
    decam_postage_stamp(img_path, binsize=120, plot_path=plot_path, save_path=None, vrange=[-vrange, vrange], dr8=False, median=True,
        blob_mask=False, ood_mask=True, show=False)

with Pool(processes=n_processes) as pool:
    res = pool.map(wrapper, np.arange(len(exp)))


