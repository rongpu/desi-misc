# cd /dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/dr11/south/tractor

import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool


def filter_tractor_cat(fn):
    cat = Table(fitsio.read(fn, columns=columns))
    # print(len(cat))

    mask = (cat['nobs_g']>=1) & (cat['nobs_r']>=1) & (cat['nobs_z']>=1)
    cat = cat[mask]

    mask = cat['brick_primary']==True
    cat = cat[mask]

    maskbits = [1, 5, 6, 7, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['maskbits'] & 2**bit)==0
    cat = cat[mask_clean]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cat['gmag'] = 22.5 - 2.5*np.log10(np.clip(cat['flux_g'], 1e-7, None))
        cat['rmag'] = 22.5 - 2.5*np.log10(np.clip(cat['flux_r'], 1e-7, None))
        cat['zmag'] = 22.5 - 2.5*np.log10(np.clip(cat['flux_z'], 1e-7, None))
        cat['gfibermag'] = 22.5 - 2.5*np.log10(np.clip(cat['fiberflux_g'], 1e-7, None))
        cat['rfibermag'] = 22.5 - 2.5*np.log10(np.clip(cat['fiberflux_r'], 1e-7, None))
        cat['zfibermag'] = 22.5 - 2.5*np.log10(np.clip(cat['fiberflux_z'], 1e-7, None))

    mask = (cat['flux_g']<=0) & (cat['rfibermag']<22) & (cat['zfibermag']<22)
    cat['outlier_g'] = mask.copy()
    mask = (cat['flux_r']<=0) & (cat['gfibermag']<22) & (cat['zfibermag']<22)
    cat['outlier_r'] = mask.copy()
    mask = (cat['flux_z']<=0) & (cat['gfibermag']<22) & (cat['rfibermag']<22)
    cat['outlier_z'] = mask.copy()

    mask = (cat['flux_g']==0) & (cat['flux_r']==0) & (cat['flux_z']==0)
    cat['outlier_zeros'] = mask.copy()

    mask = cat['outlier_g'] | cat['outlier_r'] | cat['outlier_z'] | cat['outlier_zeros']
    cat = cat[mask]

    if len(cat)==0:
        return None

    return cat


columns = ['brick_primary', 'type', 'ra', 'dec', 'flux_g', 'flux_r', 'flux_z', 'fiberflux_g', 'fiberflux_r', 'fiberflux_z', 'nobs_g', 'nobs_r', 'nobs_z', 'maskbits']

fns = glob.glob('/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/dr11/south/tractor/*/tractor-*.fits')
# fns = glob.glob('*/tractor-*.fits')
print(len(fns))

n_process = 256
with Pool(processes=n_process) as pool:
    res = pool.map(filter_tractor_cat, fns)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res)
print(len(cat))

from desiutil import brick
bricks = brick.Bricks(bricksize=0.25)
cat['brickname'] = bricks.brickname(cat["ra"], cat["dec"])
cat['brickid'] = bricks.brickid(cat["ra"], cat["dec"])

cat.write('/pscratch/sd/r/rongpu/tmp/dr11_negative_flux_outliers.fits')
