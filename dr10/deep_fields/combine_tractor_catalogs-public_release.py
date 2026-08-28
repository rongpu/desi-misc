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


fns = sorted(glob.glob('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/tractor/*/tractor-1*.fits'))

print(len(fns))


def read_catalog(index):
    fn = fns[index]
    cat = Table(fitsio.read(fn))
    cat = cat[cat['brick_primary']]

    return cat


n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(read_catalog, np.arange(len(fns)))
cat = vstack(res)

sweep = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/sweep/10.1/sweep-040p010-045p015.fits', rows=[0]))
sweep_extra = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/sweep/10.1-extra/sweep-040p010-045p015-ex.fits', rows=[0]))
sweep_lc = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr10/south/sweep/10.1-lightcurves/sweep-040p010-045p015-lc.fits', rows=[0]))

tractor_all_col = np.array(cat.colnames)

sweep_col = np.char.lower(sweep.colnames)
extra_col = np.char.lower(sweep_extra.colnames)
lc_col = np.char.lower(sweep_lc.colnames)

# Should be all empty:
print(sweep_col[~np.in1d(sweep_col, tractor_all_col)])
print(extra_col[~np.in1d(extra_col, tractor_all_col)])
print(lc_col[~np.in1d(lc_col, tractor_all_col)])

tractor = cat[list(tractor_all_col[np.in1d(tractor_all_col, sweep_col)])].copy()
tractor_extra = cat[list(tractor_all_col[np.in1d(tractor_all_col, extra_col)])].copy()
tractor_lc = cat[list(tractor_all_col[np.in1d(tractor_all_col, lc_col)])].copy()

saved_cols = np.unique(np.concatenate([tractor.colnames, tractor_extra.colnames, tractor_lc.colnames]))
assert np.all(saved_cols==np.unique(tractor_all_col))

tractor.write('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/catalogs/tractor-combined.fits')
tractor_extra.write('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/catalogs/tractor-combined-extra.fits')
tractor_lc.write('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/catalogs/tractor-combined-lc.fits')

