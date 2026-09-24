from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool


lc_columns = ['lc_flux_w1', 'lc_flux_w2', 'lc_flux_ivar_w1', 'lc_flux_ivar_w2', 'lc_nobs_w1', 'lc_nobs_w2', 'lc_fracflux_w1', 'lc_fracflux_w2', 'lc_rchisq_w1', 'lc_rchisq_w2', 'lc_mjd_w1', 'lc_mjd_w2', 'lc_epoch_index_w1', 'lc_epoch_index_w2']

fns = glob.glob('/pscratch/sd/r/rongpu/ibis-dr1/tractor/*/tractor-*.fits')
print(len(fns))

def read_catalog(index):
    fn = fns[index]
    cat = Table(fitsio.read(fn))
    cat.remove_columns(lc_columns)
    cat = cat[cat['brick_primary']]
    return cat

n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(read_catalog, np.arange(len(fns)))
cat = vstack(res)
print(len(cat))

mask = cat['ra']<100
cat_xmm = cat[mask].copy()
print('XMM', len(cat_xmm))
cat_xmm.rename_column('ls_id_dr11', 'ibis_id')
cat_xmm.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_xmm.fits')

mask = cat['ra']>100
cat_cosmos = cat[mask].copy()
print('COSMOS', len(cat_cosmos))
cat_cosmos.rename_column('ls_id_dr11', 'ibis_id')
cat_cosmos.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_cosmos.fits')
