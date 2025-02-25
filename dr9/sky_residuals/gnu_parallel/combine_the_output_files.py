# This code is too slow to finish running in 4 hours on an interactive node

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

field = 'south'
# field = 'north'

fns = glob.glob('/global/u2/r/rongpu/data/temp/survey_bricks_optical_sky_residuals-dr9_{}_*.fits'.format(field))
print(fns)

cat = []
for fn in fns:
    cat.append(Table.read(fn))

cat = vstack(cat)
cat.write('/global/cfs/cdirs/desi/users/rongpu/dr9/sky_residuals/survey_bricks_optical_sky_residuals-dr9_{}.fits'.format(field))