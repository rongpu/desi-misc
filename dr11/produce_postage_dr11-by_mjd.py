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
from decam_postage_stamps import decam_plot


n_processes = 128

img_dir = '/dvs_ro/cfs/cdirs/cosmo/staging'

# plot_dir = '/global/cfs/cdirs/cosmo/www/temp/rongpu/dr11/postage_stamps_by_plver'
plot_dir = '/pscratch/sd/r/rongpu/dr11/postage_stamps_by_mjd'

ccd = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/survey-ccds-decam-dr11-merged-trim-temp.fits'))
print(len(ccd))

mask = ccd['filter']=='z'
ccd = ccd[mask]
print(len(ccd))

mask = ccd['ccd_cuts']==0
ccd = ccd[mask]
print(len(ccd))

from astropy.coordinates import SkyCoord
c = SkyCoord(ccd['ra_bore'], ccd['dec_bore'], unit='deg').galactic
ccd['l'], ccd['b'] = c.l.to_value('deg'), c.b.to_value('deg')
mask = np.abs(ccd['b'])>20
ccd = ccd[mask]
print(len(ccd))

mask = ccd['dec_bore']>-60
ccd = ccd[mask]
print(len(ccd))

np.random.seed(88)
idx = np.random.choice(len(ccd), size=400, replace=False)
print(len(idx))

def wrapper(index):
    img_path = os.path.join(img_dir, ccd['image_filename'][index])
    band = ccd['filter'][index]
    plver = ccd['plver'][index]
    expnum = ccd['expnum'][index]
    mjd = ccd['mjd_obs'][index]
    plot_path = os.path.join(plot_dir, '{}_{}_{}.png'.format(str(int(mjd)), plver, expnum))
    # print(plot_path)
    ccdzpt, exptime = ccd['ccdzpt'][index], ccd['exptime'][index]
    vrange = np.array([-1., 1.]) * (0.5 * 0.262**2) * (10**((ccdzpt-22.5)/2.5) * exptime)  # normalize to +- 0.5 nanomaggies per arcsec^2
    decam_plot(os.path.join(img_path), cmap='gray', plot_path=plot_path, vrange=vrange)

with Pool(processes=n_processes) as pool:
    res = pool.map(wrapper, idx, chunksize=1)


