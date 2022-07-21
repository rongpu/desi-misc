from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio


field_index = 1

field_names = ['s12', 'x123', 'c123', 'e12']
radec_limits = [[38.9, 45.1, -3.3, 2.3],
                [32.2, 38.8, -8.7, -2.3],
                [50, 56.9, -31.5, -24.8],
                [4.7, 12.7, -46.3, -40.7]]

bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields-des-sn.fits'))
print(len(bricks))

ramin, ramax, decmin, decmax = radec_limits[field_index]
mask = (bricks['RA']>ramin) & (bricks['RA']<ramax) & (bricks['DEC']>decmin) & (bricks['DEC']<decmax)
bricks = bricks[mask]
print(len(bricks))

# Select bricks within 0.9 deg of the center of X3
sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord
deep_ra = np.array([36.4500])
deep_dec = np.array([-4.6000])
search_radius = 0.9*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, bricks['RA'], bricks['DEC'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)

bricks = bricks[idx2]

np.random.seed(98171)
idx = np.random.choice(len(bricks), size=len(bricks), replace=False)  # shuffle
bricks = bricks[idx]

# fn = '/global/u2/r/rongpu/temp/run/tasks-runbrick-x3.sh'
# with open(fn, 'w') as f:
#     for brickname in bricks['BRICKNAME']:
#         if not os.path.isfile('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/tractor/{}/tractor-{}.fits'.format(brickname[:3], brickname)):
#             f.write('shifter --image docker:legacysurvey/legacypipe:DR10.0.0 ./runbrick-cosmos.sh {} \n'.format(brickname))

fn = '/global/u2/r/rongpu/temp/qdo/bricks-x3.txt'
with open(fn, 'w') as f:
    for brickname in bricks['BRICKNAME']:
        if not os.path.isfile('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/tractor/{}/tractor-{}.fits'.format(brickname[:3], brickname)):
            f.write(brickname+'\n')
