# Print the list of brick names for runbrick-splinesky.sh

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

# Load CCD list
ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-defringed.fits'))
print(len(ccd))

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx].copy()
print(len(exp))

mask = exp['ra_bore']>130
exp = exp[mask]
print(len(exp))

mask = exp['filter']!='Y'
exp = exp[mask]
print(len(exp))

command_output_path = '/global/u2/r/rongpu/temp/tractor/splinesky_cosmos_commands.sh'
with open(command_output_path, 'w') as f:
    for image_filename in exp['image_filename']:
        f.write('shifter --image docker:legacysurvey/legacypipe:DR10.0.0 ./runcalibs-sky.sh {}\n'.format(image_filename))
