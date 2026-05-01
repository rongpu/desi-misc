# Combine the CCDs from ibis/deep_field_subsets/depth_cut.py

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

params = {'legend.fontsize': 'medium',
          'axes.labelsize': 'medium',
          'axes.titlesize': 'medium',
          'xtick.labelsize': 'medium',
          'ytick.labelsize': 'medium',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']
field_names = ['COSMOS', 'XMM-LSS']

ccd_stack = []
for field_index in [0, 1]:
    for band in ibis_filters:
        field = field_names[field_index]
        print(field, band)
        ccd = Table(fitsio.read(f'/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/misc/survey-ccds-ibis-dr1-subset_25.0_{field}_{band}.fits'))
        # ccd = Table(fitsio.read(f'/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/misc/survey-ccds-ibis-dr1-subset_25.0_{field}_{band}-1.fits'))
        # ccd = Table(fitsio.read(f'/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/misc/survey-ccds-ibis-dr1-subset_24.8_{field}_{band}.fits'))
        # ccd = Table(fitsio.read(f'/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/misc/survey-ccds-ibis-dr1-subset_25.2_{field}_{band}.fits'))
        print(len(ccd))
        ccd_stack.append(ccd)
ccd = vstack(ccd_stack)
print(len(ccd))
ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/survey-ccds-ibis-dr1-subset-25.0.fits')
# ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/survey-ccds-ibis-dr1-subset-25.0-1.fits')
# ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/survey-ccds-ibis-dr1-subset-24.8.fits')
# ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/survey-ccds-ibis-dr1-subset-25.2.fits')
