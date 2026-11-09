# To-do: add filler targets

import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

cat1 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-ibis_deep.fits'))
print(len(cat1))
cat1 = cat1[cat1['new_obs']]
print(len(cat1))

cat2 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-ibis_wide.fits'))
print(len(cat2))
cat2 = cat2[cat2['new_obs']]
print(len(cat2))

cat1 = cat1[['ibis_id_dr1', 'ra', 'dec']]
cat2 = cat2[['ibis_id_dr1', 'ra', 'dec']]
cat = vstack([cat1, cat2])
print(len(cat))

_, idx = np.unique(cat['ibis_id_dr1'], return_index=True)
cat = cat[idx]

mask = np.in1d(cat['ibis_id_dr1'], cat1['ibis_id_dr1'])
cat['ibis_deep_lae'] = mask.copy()
mask = np.in1d(cat['ibis_id_dr1'], cat2['ibis_id_dr1'])
cat['ibis_wide_lae'] = mask.copy()

cat['high_priority'] = cat['ibis_deep_lae'].copy()
cat['mid_priority'] = ~cat['high_priority']
assert np.all(cat['ibis_wide_lae'] & (~cat['ibis_deep_lae'])==cat['mid_priority'])

print(np.sum(cat['high_priority']), np.sum(cat['high_priority'])/len(cat))
print(np.sum(cat['mid_priority']), np.sum(cat['mid_priority'])/len(cat))

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-combined.fits')
