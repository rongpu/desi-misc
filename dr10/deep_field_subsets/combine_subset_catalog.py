from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

# Columns to discard
apflux_columns = ['apflux_g', 'apflux_r', 'apflux_z', 'apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z', 'apflux_ivar_g', 'apflux_ivar_r', 'apflux_ivar_z', 'apflux_masked_g', 'apflux_masked_r', 'apflux_masked_z', 'apflux_w1', 'apflux_w2', 'apflux_w3', 'apflux_w4', 'apflux_resid_w1', 'apflux_resid_w2', 'apflux_resid_w3', 'apflux_resid_w4', 'apflux_ivar_w1', 'apflux_ivar_w2', 'apflux_ivar_w3', 'apflux_ivar_w4']
lc_columns = ['lc_flux_w1', 'lc_flux_w2', 'lc_flux_ivar_w1', 'lc_flux_ivar_w2', 'lc_nobs_w1', 'lc_nobs_w2', 'lc_fracflux_w1', 'lc_fracflux_w2', 'lc_rchisq_w1', 'lc_rchisq_w2', 'lc_mjd_w1', 'lc_mjd_w2', 'lc_epoch_index_w1', 'lc_epoch_index_w2']
apflux_columns_dr10 = ['apflux_g', 'apflux_r', 'apflux_i', 'apflux_z', 'apflux_resid_g', 'apflux_resid_r', 'apflux_resid_i', 'apflux_resid_z', 'apflux_ivar_g', 'apflux_ivar_r', 'apflux_ivar_i', 'apflux_ivar_z', 'apflux_masked_g', 'apflux_masked_r', 'apflux_masked_i', 'apflux_masked_z', 'apflux_w1', 'apflux_w2', 'apflux_w3', 'apflux_w4', 'apflux_resid_w1', 'apflux_resid_w2', 'apflux_resid_w3', 'apflux_resid_w4', 'apflux_ivar_w1', 'apflux_ivar_w2', 'apflux_ivar_w3', 'apflux_ivar_w4']

# Rongpu's COSMOS subsets
fns = glob.glob('/global/cfs/cdirs/desi/users/rongpu/data/deep_field_subsets/cosmos/*/tractor/*/tractor-*.fits')
cat = []
for fn in fns:
    tmp = Table(fitsio.read(fn))
    ii = fn.find('/cosmos/')
    jj = fn.find('/tractor/')
    tmp['sub'] = fn[ii+8:jj]
    tmp.remove_columns(lc_columns)
    tmp.remove_columns(apflux_columns)
    tmp = tmp[tmp['brick_primary']]
    cat.append(tmp)
cat = vstack(cat)
cat.write('/global/cfs/cdirs/desi/users/rongpu/data/deep_field_subsets/catalogs/cosmos_subsets_rongpu_dr10.fits')

# Dustin's DR9 COSMOS subsets
fns = glob.glob('/global/cfs/cdirs/desi/users/rongpu/dr9/dr9-cosmos-subs/*/tractor/*/*.fits')
cat = []
for fn in fns:
    tmp = Table(fitsio.read(fn))
    ii = fn.find('/tractor/')
    tmp['sub'] = fn[ii-2:ii]
    tmp.remove_columns(lc_columns)
    tmp.remove_columns(apflux_columns)
    tmp = tmp[tmp['brick_primary']]
    cat.append(tmp)
cat = vstack(cat)
cat.write('/global/cfs/cdirs/desi/users/rongpu/data/deep_field_subsets/catalogs/cosmos_subsets_dstn_dr9.fits')

# Dustin's DR10 COSMOS subsets
fns = glob.glob('/global/cscratch1/sd/dstn/dr10-cosmos-subs/*/tractor/*/*.fits')
cat = []
for fn in fns:
    tmp = Table(fitsio.read(fn))
    ii = fn.find('/tractor/')
    tmp['sub'] = fn[ii-2:ii]
    tmp.remove_columns(lc_columns)
    tmp.remove_columns(apflux_columns_dr10)
    tmp = tmp[tmp['brick_primary']]
    cat.append(tmp)
cat = vstack(cat)
cat.write('/global/cfs/cdirs/desi/users/rongpu/data/deep_field_subsets/catalogs/cosmos_subsets_dstn_dr10.fits')
