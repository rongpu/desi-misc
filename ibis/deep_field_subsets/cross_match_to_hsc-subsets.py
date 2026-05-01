from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


# subsets = ['subset-25.0', 'subset-25.0-1', 'subset-24.8', 'subset-25.2']
subsets = ['subset-25.0', 'subset-25.0-1', 'subset-24.8']

############################################ HSC-wide ############################################
hsc_fns = ['/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-wide-xmm-lss-no_mag_limit-reduced.fits',
       '/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-wide-cosmos-no_mag_limit-reduced.fits']

for hsc_fn in hsc_fns:
    hsc_wide = Table(fitsio.read(hsc_fn))
    print(len(hsc_wide))

    for subset in subsets:

        if 'xmm-lss' in hsc_fn:
            ibis_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/{}/tractor_xmm.fits'.format(subset)
        else:
            ibis_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/{}/tractor_cosmos.fits'.format(subset)
        cat = Table(fitsio.read(ibis_fn, columns=['ra', 'dec', 'ibis_id']))
        print(len(cat))
        cat['id'] = np.arange(len(cat))

        idx1, idx2, d2d, d_ra, d_dec = match_coord(hsc_wide['ra'], hsc_wide['dec'], cat['ra'], cat['dec'], search_radius=1., plot_q=False)
        print(len(idx1)/len(cat))
        print(len(idx1)/len(hsc_wide))
        hsc = hsc_wide[idx1].copy()
        hsc['id'] = idx2

        tt = cat[['id', 'ibis_id']].copy()
        tt = join(tt, hsc, keys='id', join_type='left').filled(-99)
        tt.remove_column('id')
        tt.write(ibis_fn.replace('.fits', '-hsc_wide.fits'))


############################################ HSC-deep ############################################

hsc_deep = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-dud-rev-no_mag_limit-reduced.fits'))
print(len(hsc_deep))

for subset in subsets:

    for ibis_fn in ['/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/{}/tractor_xmm.fits'.format(subset),
                    '/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/{}/tractor_cosmos.fits'.format(subset)]:
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
