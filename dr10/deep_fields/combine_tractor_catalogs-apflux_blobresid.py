# Combine the tractor FITS catalogs
# Discard the light curve and aperture flux columns

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool

# output_path = '/global/cfs/cdirs/desi/users/rongpu/data/decam_deep_fields/cosmos.fits'
output_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/catalogs/cosmos_apflux_blobresid.fits'
fns = glob.glob('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/tractor/*/tractor-1*.fits')
# fns = glob.glob('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/tractor/*/tractor-1*.fits')

# output_path = '/global/cfs/cdirs/desi/users/rongpu/data/decam_deep_fields/des_sn.fits'
# fns = glob.glob('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/tractor/*/tractor-0*.fits')

print(len(fns))

columns = ['brick_primary', 'apflux_blobresid_g', 'apflux_blobresid_r', 'apflux_blobresid_i', 'apflux_blobresid_z']


def read_catalog(index):
    fn = fns[index]
    cat = Table(fitsio.read(fn, columns=columns))
    cat['n_sources_all'] = len(cat)
    cat = cat[cat['brick_primary']]

    cat.remove_column('brick_primary')

    return cat


n_processes = 256
with Pool(processes=n_processes) as pool:
    res = pool.map(read_catalog, np.arange(len(fns)))
cat = vstack(res)

cat.write(output_path, overwrite=True)
