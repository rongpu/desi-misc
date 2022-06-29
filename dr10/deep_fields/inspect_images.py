from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/useful/'))
from decam_postage_stamps import decam_plot


n_processes = 128

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1.fits'))
print(len(cat))

mask = cat['ra']>130
cat = cat[mask]

deep_ra = np.array([54.2743, 54.2743, 52.6484, 34.4757, 35.6645, 36.4500, 42.8200, 41.1944, 7.8744, 9.5000, 150.1166])
deep_dec = np.array([-27.1116, -29.0884, -28.1000, -4.9295, -6.4121, -4.6000, 0.0000, -0.9884, -43.0096, -43.9980, 2.2058])

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

search_radius = 1.3*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, cat['ra'], cat['dec'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
print(len(np.unique(cat['expnum'][idx2])))

cat = cat[idx2]

_, idx = np.unique(cat['expnum'], return_index=True)
exp = cat[idx]
print(len(exp), len(cat)/len(exp))


def make_plot(index):
    fn = '/global/cfs/cdirs/cosmo/staging/' + exp['image_filename'][index]
    expnum = exp['expnum'][index]
    band = exp['filter'][index]
    # plot_path = '/pscratch/sd/r/rongpu/dr10dev/deep_field/postage_stamps/cosmos/{}_{}.png'.format(expnum, band)
    # decam_plot(fn, plot_path=plot_path)
    plot_path = '/pscratch/sd/r/rongpu/dr10dev/deep_field/postage_stamps/cosmos/{}_{}_ood.png'.format(expnum, band)
    decam_plot(fn, plot_path=plot_path, ood_mask=True)
    # fn = '/global/cfs/cdirs/cosmo/staging/' + exp['image_filename'][index].replace('_ooi_', '_ood_')
    # decam_plot(fn, cmap='gray_r', vrange=[0, 1], median=False, subtract_median_sky=False)


print('Start!')
time_start = time.time()

with Pool(processes=n_processes) as pool:
    res = pool.map(make_plot, np.arange(len(exp)))

print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))
