# Print the list of brick names for runbrick-splinesky.sh

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

bricks_path = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-bricks.fits.gz'
bricks = Table(fitsio.read(bricks_path, columns=['RA', 'DEC']))

idx = np.arange(len(bricks))
mask = (bricks['DEC']>-90) & (bricks['DEC']<90)
mask &= bricks['RA']>130  # COSMOS
bricks = bricks[mask]
idx = idx[mask]

search_radius = 3.0*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, bricks['RA'], bricks['DEC'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
idx2 = np.sort(idx2)
idx = idx[idx2]

bricks = Table(fitsio.read(bricks_path, rows=idx))

command_output_path = '/global/u2/r/rongpu/temp/tractor/splinesky_cosmos_commands.sh'
with open(command_output_path, 'w') as f:
    for brickname in bricks['BRICKNAME']:
        f.write('./runbrick-splinesky.sh {}\n'.format(brickname))
