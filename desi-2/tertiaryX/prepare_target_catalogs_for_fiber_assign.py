from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

fns = [
'/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/tertiary54_lae_targets.fits',
'/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/tertiary54_lae_filler_targets.fits',
'/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/LSST_Y4_CC_RF_COSMOS_lbg_targets_unique.fits',
]

for fn in fns:

    cat = Table(fitsio.read(fn))

    for col in ['ra', 'dec']:
        if col in cat.colnames:
            cat.rename_column(col, col.upper())

    cat['PMRA'] = np.zeros(cat['RA'].shape)
    cat['PMDEC'] = np.zeros(cat['RA'].shape)
    cat['REF_EPOCH'] = np.ones(cat['RA'].shape)* 2015.5
    cat = cat[['RA', 'DEC', 'PMRA', 'PMDEC', 'REF_EPOCH']]

    cat.write(fn.replace('/tertiary54/', '/tertiary54/radec_only/').replace('.fits', '_radec.fits'))
