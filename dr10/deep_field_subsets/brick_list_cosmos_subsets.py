from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


deep_ra = np.array([150.1166])
deep_dec = np.array([2.2058])

bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-bricks.fits.gz'))
print(len(bricks))

search_radius = 0.95*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, bricks['RA'], bricks['DEC'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
idx2 = np.sort(idx2)
bricks = bricks[idx2]

np.random.seed(98175)
idx = np.random.choice(len(bricks), size=len(bricks), replace=False)  # shuffle
bricks = bricks[idx]

fn = '/global/u2/r/rongpu/temp/subs/bricks-cosmos-subs.txt'
# with open(fn, 'w') as f:
#     for brickname in bricks['BRICKNAME']:
#         if not os.path.isfile('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/tractor/{}/tractor-{}.fits'.format(brickname[:3], brickname)):
#             f.write('shifter --image docker:legacysurvey/legacypipe:DR10.0.0 ./runbrick-cosmos.sh {} \n'.format(brickname))
with open(fn, 'w') as f:
    for brickname in bricks['BRICKNAME']:
        f.write('{}\n'.format(brickname))
