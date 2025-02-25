# Replace the skyscales for the newly reprocessed V4.9 images

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits

skyscale1 = Table.read('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds.fits')
skyscale2 = Table.read('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_cp_v4.9.fits')
print(len(skyscale1))
print(len(skyscale2))

skyscale1['ccdname'] = np.char.strip(skyscale1['ccdname']) # strip spaces
skyscale2['ccdname'] = np.char.strip(skyscale2['ccdname']) # strip spaces
skyscale1.remove_column('image_hdu')

mask = np.in1d(skyscale1['expnum'], skyscale2['expnum'])
skyscale1 = skyscale1[~mask]
print(len(skyscale1))

skyscale_new = vstack([skyscale1, skyscale2])

# Sanity check: all the CCDs should be unique
skyscale_new['ccd_id'] = np.char.add(np.array(skyscale_new['expnum'], dtype=str), np.array(skyscale_new['ccdname'], dtype=str))
if len(np.unique(skyscale_new['ccd_id']))!=len(skyscale_new):
    raise ValueError

skyscale_new.sort('ccd_id')
skyscale_new.remove_column('ccd_id')

skyscale_new.write('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_new.fits')
