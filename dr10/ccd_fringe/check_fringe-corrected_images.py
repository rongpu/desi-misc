from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool


n_processes = 128

image_output_dir = '/pscratch/sd/r/rongpu/dr10dev/fringe_corrected_images'

# Load CCD list
ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1.fits'))
print(len(ccd))

mask = ccd['filter']=='z'  # include only z-band images
ccd = ccd[mask]

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx].copy()
print(len(exp))

# Shuffle
np.random.seed(781234)
idx = np.random.choice(len(exp), size=len(exp), replace=False)
exp = exp[idx]


def check_corrected_image(index):
    fn = os.path.join(image_output_dir, exp['image_filename'][index])
    if os.path.isfile(fn):
        n_images = len(fitsio.FITS(fn)) - 1
        if n_images!=exp['n_ccds'][index]:
            print('Missing CCDs: {}/{}  {} {}'.format(n_images, exp['n_ccds'][index], exp['expnum'][index], exp['image_filename'][index]))
    else:
        print('Missing file: {} {}'.format(exp['expnum'][index], exp['image_filename'][index]))


print('Start!')
time_start = time.time()

with Pool(processes=n_processes) as pool:
    res = pool.map(check_corrected_image, np.arange(len(exp)))

print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))
