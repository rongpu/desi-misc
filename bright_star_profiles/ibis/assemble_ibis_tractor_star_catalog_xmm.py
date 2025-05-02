from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

fns = glob.glob('/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis3/tractor/*/tractor-*.fits')
print(len(fns))
fns += glob.glob('/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis3-bail/tractor/*/tractor-*.fits')
print(len(fns))
fns.sort()
print(len(fns))

cat_stack = []
for fn in fns:
    cat = Table(fitsio.read(fn))
    # print(len(cat))
    mask = cat['gaia_phot_g_mean_mag']!=0
    mask &= cat['brick_primary']==True
    mask &= cat['maskbits'] & 2**10 == 0  # BAILOUT
    mask &= (cat['gaia_phot_bp_mean_mag']!=0) & (cat['gaia_phot_rp_mean_mag']!=0)
    cat = cat[mask]
    # print(len(cat))
    
    if 'ibis3-bail' in fn:
        cat['bailout_brick'] = True
    else:
        cat['bailout_brick'] = False

    cat_stack.append(cat)
cat = vstack(cat_stack)
print(len(cat))

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    cat['m411'] = 22.5 - 2.5*np.log10(cat['flux_M411'])
    cat['m438'] = 22.5 - 2.5*np.log10(cat['flux_M438'])
    cat['m464'] = 22.5 - 2.5*np.log10(cat['flux_M464'])
    cat['m490'] = 22.5 - 2.5*np.log10(cat['flux_M490'])
    cat['m517'] = 22.5 - 2.5*np.log10(cat['flux_M517'])

mask_bad = np.full(len(cat), False)
for band in ['m411', 'm438', 'm464', 'm490', 'm517']:
    mask_bad |= cat['flux_ivar_'+band.upper()]==0
print(np.sum(mask_bad))

cat = cat[~mask_bad]
print(len(cat))

cat.write('/pscratch/sd/r/rongpu/ibis/tractor_xmm_stars.fits')
