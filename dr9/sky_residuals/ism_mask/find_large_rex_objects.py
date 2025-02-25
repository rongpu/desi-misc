# Find large objects that are mostly spurious sources from emission clouds

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool


def find_large_objects(sweep_fn):

    columns = ['RELEASE', 'BRICKID', 'OBJID', 'RA', 'DEC', 'TYPE', 'SHAPE_R', 'MASKBITS']
    cat = Table(fitsio.read(sweep_fn, columns=columns))

    # ext_columns = ['APFLUX_BLOBRESID_G', 'APFLUX_BLOBRESID_R', 'APFLUX_BLOBRESID_Z']
    # sweep_ext_fn = sweep_fn.replace('/9.0/', '/9.0-extra/').replace('.fits', '-ex.fits')
    # cat_ext = Table(fitsio.read(sweep_ext_fn, columns=ext_columns))
    # cat = hstack([cat, cat_ext])

    mask = cat['SHAPE_R']>2
    # mask &= cat['MASKBITS']==0
    mask &= cat['TYPE']=='REX'
    print(os.path.basename(sweep_fn), np.sum(mask))

    if np.sum(mask)==0:
        return None
    else:
        cat = cat[mask].copy()
        if 'north' in sweep_fn:
            cat['PHOTSYS'] = 'N'
        else:
            cat['PHOTSYS'] = 'S'
        return cat


sweep_fns = glob.glob('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/*/sweep/9.0/*.fits')
sweep_fns.sort()

n_process = 64
with Pool(processes=n_process) as pool:
    res = pool.map(find_large_objects, sweep_fns)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res)
cat.write('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/large_rex_objects.fits', overwrite=True)
