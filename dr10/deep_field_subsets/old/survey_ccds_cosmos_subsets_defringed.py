# Switch to the defringed images for z band

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

suffixes = ['dr9-ccds'] + list(range(11))

for suffix in suffixes:

    ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/survey-ccds-dr10-v6-subset-{}.fits'.format(suffix)))
    max_str_length = np.max([len(ccd['image_filename'][index]) for index in range(len(ccd))])

    if max_str_length+10>120:
        raise ValueError('The new image_filename exceeds the current string length.')

    mask = ccd['filter']=='z'
    ccd['image_filename'][mask] = np.char.replace(ccd['image_filename'][mask], 'decam/CP', 'decam-defringed/CP')

    ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/defringed/survey-ccds-dr10-v6-subset-{}-defringed.fits'.format(suffix), overwrite=True)
