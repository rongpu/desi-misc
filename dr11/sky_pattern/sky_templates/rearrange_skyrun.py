# rearrange DR11 skyrun run numbers so that it can be combined with DR9 templates

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

skyrun = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/schlafly/decals/skyrunsdr11-v3.fits'))
print('DR11', len(skyrun))

mask = skyrun['filter']=='z'
skyrun = skyrun[mask]
print('DR11', len(skyrun))


dr9 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr9/sky_pattern/skyrunsgoodcountexpnumv48dr8.fits'))
print('DR9', len(dr9))

mask = dr9['filter']=='z'
dr9 = dr9[mask]
print('DR9', len(dr9))

mask = skyrun['mjd_obs']>dr9['mjd_obs'].max()
skyrun = skyrun[mask]
print('DR11', len(skyrun))

skyrun['run'] = skyrun['run'] - skyrun['run'].min() + 1 + dr9['run'].max()

skyrun.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/skyrunsdr11-v3-rearrange.fits')

