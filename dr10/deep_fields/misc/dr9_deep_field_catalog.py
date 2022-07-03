from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

cat_stack = []
lc_columns = ['lc_flux_w1', 'lc_flux_w2', 'lc_flux_ivar_w1', 'lc_flux_ivar_w2', 'lc_nobs_w1', 'lc_nobs_w2', 'lc_fracflux_w1', 'lc_fracflux_w2', 'lc_rchisq_w1', 'lc_rchisq_w2', 'lc_mjd_w1', 'lc_mjd_w2', 'lc_epoch_index_w1', 'lc_epoch_index_w2']
apflux_columns = ['apflux_g', 'apflux_r', 'apflux_z', 'apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z', 'apflux_blobresid_g', 'apflux_blobresid_r', 'apflux_blobresid_z', 'apflux_ivar_g', 'apflux_ivar_r', 'apflux_ivar_z', 'apflux_masked_g', 'apflux_masked_r', 'apflux_masked_z', 'apflux_w1', 'apflux_w2', 'apflux_w3', 'apflux_w4', 'apflux_resid_w1', 'apflux_resid_w2', 'apflux_resid_w3', 'apflux_resid_w4', 'apflux_ivar_w1', 'apflux_ivar_w2', 'apflux_ivar_w3', 'apflux_ivar_w4']
fns = glob.glob('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9.1.1/tractor/*/*.fits')
for fn in fns:
    cat = Table(fitsio.read(fn))
    cat.remove_columns(lc_columns)
    cat.remove_columns(apflux_columns)
    mask = cat['brick_primary'].copy()
    cat = cat[mask]
    cat_stack.append(cat)
cat_stack = vstack(cat_stack)
cat_stack.write('/global/cfs/cdirs/desi/users/rongpu/dr9/misc/cosmos_deep_dield_9.1.1.fits')
