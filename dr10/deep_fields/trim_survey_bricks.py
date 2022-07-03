# Trim survey-bricks-10x10.fits and survey-bricks-4x4.fits so their RA, DEC are within 3 degrees of the field centers

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


deep_ra = np.array([54.2743, 54.2743, 52.6484, 34.4757, 35.6645, 36.4500, 42.8200, 41.1944, 7.8744, 9.5000, 150.1166])
deep_dec = np.array([-27.1116, -29.0884, -28.1000, -4.9295, -6.4121, -4.6000, 0.0000, -0.9884, -43.0096, -43.9980, 2.2058])

data_dir = '/global/cfs/cdirs/cosmo/work/users/rongpu'

for fn in ['survey-bricks-4x4.fits', 'survey-bricks-10x10.fits']:
    path = os.path.join(data_dir, fn)
    print(path)
    bricks = Table(fitsio.read(path, columns=['RA', 'DEC']))

    idx = np.arange(len(bricks))
    mask = (bricks['DEC']>-90) & (bricks['DEC']<90)
    bricks = bricks[mask]
    idx = idx[mask]

    search_radius = 3.0*3600.  # arcsec
    idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, bricks['RA'], bricks['DEC'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
    idx2 = np.sort(idx2)
    idx = idx[idx2]

    bricks = Table(fitsio.read(path, rows=idx))
    bricks.write(path.replace('.fits', '-decam-deep-fields.fits'), overwrite=True)
