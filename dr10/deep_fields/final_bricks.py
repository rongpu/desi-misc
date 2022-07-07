# Apply NCCD cut and restrict to 2.05 radius

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

from multiprocessing import Pool


deep_ra = np.array([54.2743, 54.2743, 52.6484, 34.4757, 35.6645, 36.4500, 42.8200, 41.1944, 7.8744, 9.5000, 150.1166])
deep_dec = np.array([-27.1116, -29.0884, -28.1000, -4.9295, -6.4121, -4.6000, 0.0000, -0.9884, -43.0096, -43.9980, 2.2058])

bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields.fits'))

search_radius = 2.05*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, bricks['RA'], bricks['DEC'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
bricks = bricks[idx2]
bricks.sort('BRICKID')

mask = bricks['RA']>130
bricks1 = bricks[mask].copy()
bricks1.write('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields-cosmos.fits')

mask = bricks['RA']<=130
bricks1 = bricks[mask].copy()
bricks1.write('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields-des-sn.fits')

