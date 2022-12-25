from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields-cosmos.fits'))
print(len(bricks))

np.random.seed(98175)
idx = np.random.choice(len(bricks), size=len(bricks), replace=False)  # shuffle
bricks = bricks[idx]

fn = '/global/u2/r/rongpu/temp/run/tasks-runbrick-cosmos.sh'
with open(fn, 'w') as f:
    for brickname in bricks['BRICKNAME']:
        if not os.path.isfile('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/tractor/{}/tractor-{}.fits'.format(brickname[:3], brickname)):
            f.write('shifter --image docker:legacysurvey/legacypipe:DR10.0.0 ./runbrick-cosmos.sh {} \n'.format(brickname))

# # One task per file
# counter = 0
# for brickname in bricks['BRICKNAME']:
#     if not os.path.isfile('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/tractor/{}/tractor-{}.fits'.format(brickname[:3], brickname)):
#         fn = '/global/u2/r/rongpu/temp/run/single_task/tasks-runbrick-cosmos-{}.sh'.format(counter)
#         with open(fn, 'w') as f:
#             print(fn)
#             f.write('shifter --image docker:legacysurvey/legacypipe:DR10.0.0 ./runbrick-cosmos.sh {} \n'.format(brickname))
#             counter += 1
