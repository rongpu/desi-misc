from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool


# Columns to discard
lc_columns = ['lc_flux_w1', 'lc_flux_w2', 'lc_flux_ivar_w1', 'lc_flux_ivar_w2', 'lc_nobs_w1', 'lc_nobs_w2', 'lc_fracflux_w1', 'lc_fracflux_w2', 'lc_rchisq_w1', 'lc_rchisq_w2', 'lc_mjd_w1', 'lc_mjd_w2', 'lc_epoch_index_w1', 'lc_epoch_index_w2']
apflux_columns = ['apflux_g', 'apflux_r', 'apflux_i', 'apflux_z', 'apflux_resid_g', 'apflux_resid_r', 'apflux_resid_i', 'apflux_resid_z', 'apflux_blobresid_g', 'apflux_blobresid_r', 'apflux_blobresid_i', 'apflux_blobresid_z', 'apflux_ivar_g', 'apflux_ivar_r', 'apflux_ivar_i', 'apflux_ivar_z', 'apflux_masked_g', 'apflux_masked_r', 'apflux_masked_i', 'apflux_masked_z', 'apflux_w1', 'apflux_w2', 'apflux_w3', 'apflux_w4', 'apflux_resid_w1', 'apflux_resid_w2', 'apflux_resid_w3', 'apflux_resid_w4', 'apflux_ivar_w1', 'apflux_ivar_w2', 'apflux_ivar_w3', 'apflux_ivar_w4']

# output_path = '/global/cfs/cdirs/desi/users/rongpu/data/decam_deep_fields/cosmos.fits'
output_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/catalogs/cosmos.fits'

# fns = glob.glob('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/tractor/*/tractor-1*.fits')
fns = glob.glob('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/tractor/*/tractor-1*.fits')

# output_path = '/global/cfs/cdirs/desi/users/rongpu/data/decam_deep_fields/des_sn.fits'
# fns = glob.glob('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/tractor/*/tractor-0*.fits')

print(len(fns))


def read_catalog(index):
    fn = fns[index]
    cat = Table(fitsio.read(fn))
    cat['n_sources_all'] = len(cat)
    cat = cat[cat['brick_primary']]

    cat.remove_columns(lc_columns)
    cat.remove_columns(apflux_columns)

    return cat


n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(read_catalog, np.arange(len(fns)))
cat = vstack(res)

cat.write(output_path, overwrite=True)
