from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool


n_processes = 128

# Load CCD list
ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-defringed.fits'))
print(len(ccd))

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx].copy()
print(len(exp))


def check_splinesky(index):
    fn = '/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/calib/sky/' + exp['image_filename'][index].replace('.fits.fz', '-splinesky.fits')
    if os.path.isfile(fn):
        cat = Table(fitsio.read(fn))
        if len(cat)!=exp['n_ccds'][index]:
            print('Missing CCDs: {}/{}  {} {}'.format(len(cat), exp['n_ccds'][index], exp['expnum'][index], exp['image_filename'][index]))
    else:
        print('Missing file: {} {}'.format(exp['expnum'][index], exp['image_filename'][index]))


print('Start!')
time_start = time.time()

with Pool(processes=n_processes) as pool:
    res = pool.map(check_splinesky, np.arange(len(exp)))

print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))
