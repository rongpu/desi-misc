# Create placeholder target catalogs for COSMOS fiber assignment tests

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

params = {'legend.fontsize': 'medium',
          'axes.labelsize': 'medium',
          'axes.titlesize': 'medium',
          'xtick.labelsize': 'medium',
          'ytick.labelsize': 'medium',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

ra_center, dec_center = 150.1, 2.182
ramin, ramax, decmin, decmax = ra_center-2, ra_center+2, dec_center-2, dec_center+2

lae_filler = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/xmm_highz/tertiary49_xmm_lae_targets-in_cosmos.fits'))
print(len(lae_filler))

mask = np.full(len(lae_filler), True)
plt.figure(figsize=(12, 12))
plt.plot(lae_filler['RA'][mask], lae_filler['DEC'][mask], '.', ms=0.7, alpha=0.5)
plt.grid(alpha=0.5)
plt.axis([ramax, ramin, decmin, decmax])
# plt.gca().invert_xaxis()
plt.show()

lae = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/ebina/hiz/data/tertiary49/sel_trial/cats/sel_v3.1_ibis4.fits'))
print(len(lae))

mask = np.full(len(lae), True)
plt.figure(figsize=(12, 12))
plt.plot(lae['RA'][mask], lae['DEC'][mask], '.', ms=2, alpha=1)
plt.grid(alpha=0.5)
plt.axis([ramax, ramin, decmin, decmax])
# plt.gca().invert_xaxis()
plt.show()

# Remove duplicates
idx1, idx2, d2d, d_ra, d_dec = match_coord(lae_filler['RA'], lae_filler['DEC'], lae['RA'], lae['DEC'], search_radius=1., plot_q=True)

mask = np.full(len(lae_filler), True)
mask[idx1] = False
lae_filler = lae_filler[mask]
print(len(lae_filler), np.sum(mask)/len(mask))

tmp1 = lae[['RA', 'DEC']]
tmp1['TARGET'] = 'LAE'

tmp2 = lae_filler[['RA', 'DEC']]
tmp2['TARGET'] = 'LAE_FILLER'

cat = vstack([tmp1, tmp2])
print(len(cat))

lbg1 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/misc/LSST_Y4_CC_COSMOS_lbg_targets.fits'))
lbg2 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/misc/LSST_Y4_RF_COSMOS_lbg_targets.fits'))

mask = np.full(len(lbg1), True)
plt.figure(figsize=(12, 12))
plt.plot(lbg1['RA'][mask], lbg1['DEC'][mask], '.', ms=2, alpha=1)
plt.grid(alpha=0.5)
plt.axis([ramax, ramin, decmin, decmax])
# plt.gca().invert_xaxis()
plt.show()

mask = np.full(len(lbg2), True)
plt.figure(figsize=(12, 12))
plt.plot(lbg2['RA'][mask], lbg2['DEC'][mask], '.', ms=2, alpha=1)
plt.grid(alpha=0.5)
plt.axis([ramax, ramin, decmin, decmax])
# plt.gca().invert_xaxis()
plt.show()

# Remove duplicates
idx1, idx2, d2d, d_ra, d_dec = match_coord(lbg1['RA'], lbg1['DEC'], lbg2['RA'], lbg2['DEC'], search_radius=1., plot_q=True)

mask = np.full(len(lbg1), True)
mask[idx1] = False
lbg1 = lbg1[mask]
print(len(lbg1), np.sum(mask)/len(mask))

tmp3 = lbg1[['RA', 'DEC']]
tmp3['TARGET'] = 'LBG'

tmp4 = lbg1[['RA', 'DEC']]
tmp4['TARGET'] = 'LBG'

cat = vstack([cat, tmp3, tmp4])
print(len(cat))

np.unique(cat['TARGET'], return_counts=True)

plt.figure(figsize=(12, 12))
for target in np.unique(cat['TARGET']):
    mask = cat['TARGET']==target
    if target=='LAE_FILLER':
        plt.plot(cat['RA'][mask], cat['DEC'][mask], '.', color='0.7', ms=0.1, alpha=1, label=target)
    else:
        plt.plot(cat['RA'][mask], cat['DEC'][mask], '.', ms=2, alpha=1, label=target)
plt.grid(alpha=0.5)
plt.axis([ramax, ramin, decmin, decmax])
# plt.gca().invert_xaxis()
plt.legend(markerscale=5)
plt.show()

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/cosmos_lae_lbg_placeholder.fits')

