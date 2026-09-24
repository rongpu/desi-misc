from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


############################################ HSC-wide ############################################
hsc_fns = ['/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-wide-xmm-lss-no_mag_limit-reduced.fits',
       '/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-wide-cosmos-no_mag_limit-reduced.fits']

for hsc_fn in hsc_fns:
    hsc = Table(fitsio.read(hsc_fn))
    print(len(hsc))

    if 'xmm-lss' in hsc_fn:
        ibis_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_xmm.fits'
    else:
        ibis_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_cosmos.fits'
    cat = Table(fitsio.read(ibis_fn, columns=['ra', 'dec', 'ibis_id']))
    print(len(cat))
    cat['id'] = np.arange(len(cat))

    idx1, idx2, d2d, d_ra, d_dec = match_coord(hsc['ra'], hsc['dec'], cat['ra'], cat['dec'], search_radius=1., plot_q=False)
    print(len(idx1)/len(cat))
    print(len(idx1)/len(hsc))
    hsc = hsc[idx1]
    hsc['id'] = idx2

    tt = cat[['id', 'ibis_id']].copy()
    tt = join(tt, hsc, keys='id', join_type='left').filled(-99)
    tt.remove_column('id')
    tt.write(ibis_fn.replace('.fits', '-hsc_wide.fits'))


############################################ HSC-deep ############################################

hsc_deep = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-dud-rev-no_mag_limit-reduced.fits'))
print(len(hsc_deep))

for ibis_fn in ['/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_xmm.fits',
                '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_cosmos.fits']:
    print(ibis_fn)

    cat = Table(fitsio.read(ibis_fn, columns=['ra', 'dec', 'ibis_id']))
    print(len(cat))
    cat['id'] = np.arange(len(cat))

    ramin, ramax, decmin, decmax = cat['ra'].min()-0.01, cat['ra'].max()+0.01, cat['dec'].min()-0.01, cat['dec'].max()+0.01
    mask = (hsc_deep['ra']>ramin) & (hsc_deep['ra']<ramax) & (hsc_deep['dec']>decmin) & (hsc_deep['dec']<decmax)
    hsc = hsc_deep[mask].copy()
    print(len(hsc))

    idx1, idx2, d2d, d_ra, d_dec = match_coord(hsc['ra'], hsc['dec'], cat['ra'], cat['dec'], search_radius=1., plot_q=False)
    print(len(idx1)/len(cat))
    print(len(idx1)/len(hsc))
    hsc = hsc[idx1]
    hsc['id'] = idx2

    tt = cat[['id', 'ibis_id']].copy()
    tt = join(tt, hsc, keys='id', join_type='left').filled(-99)
    tt.remove_column('id')
    tt.write(ibis_fn.replace('.fits', '-hsc_deep.fits'))
