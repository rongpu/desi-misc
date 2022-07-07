# Only unpatched PSFEx models are available for Y band
# legacypipe searches $LEGACY_SURVEY_DIR/calib/ for PSFEx files

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1.fits'))
mask = ccd['filter']=='Y'
ccd = ccd[mask]

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx]
print(len(exp))

psfex_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/calib/unpatched-psfex'
symlink_dir = '/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/runcalibs-sky/calib'

for index in range(len(exp)):
    image_filename = exp['image_filename'][index]
    psfex_filename = image_filename.replace('.fits.fz', '-psfex.fits')

    psfex_path = os.path.join(psfex_dir, psfex_filename)
    symlink_path = os.path.join(symlink_dir, psfex_filename)
    if not os.path.exists(os.path.dirname(symlink_path)):
        try:
            os.makedirs(os.path.dirname(symlink_path))
        except:
            pass
    os.symlink(psfex_path, symlink_path)
    print(symlink_path)

