# Switch to the defringed images for z band

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1.fits'))
mask = ccd['filter']=='z'
ccd = ccd[mask]

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx]
print(len(exp))

image_dir = '/global/cfs/cdirs/cosmo/staging'
new_image_dir = '/pscratch/sd/r/rongpu/dr10dev/fringe_corrected_images'
symlink_dir = '/pscratch/sd/r/rongpu/dr10dev/fringe_corrected_images_with_symlinks'

for index in range(len(exp)):
    image_filename = exp['image_filename'][index]
    if not os.path.exists(os.path.dirname(os.path.join(symlink_dir, image_filename))):
        try:
            os.makedirs(os.path.dirname(os.path.join(symlink_dir, image_filename)))
        except:
            pass
    ood_filename = image_filename.replace('_ooi_', '_ood_')
    oow_filename = image_filename.replace('_ooi_', '_oow_')
    os.symlink(os.path.join(new_image_dir, image_filename), os.path.join(symlink_dir, image_filename))
    os.symlink(os.path.join(image_dir, ood_filename), os.path.join(symlink_dir, ood_filename))
    os.symlink(os.path.join(image_dir, oow_filename), os.path.join(symlink_dir, oow_filename))

