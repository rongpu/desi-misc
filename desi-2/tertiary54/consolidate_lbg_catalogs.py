from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

# params = {'legend.fontsize': 'large',
#          'axes.labelsize': 'large',
#          'axes.titlesize':'large',
#          'xtick.labelsize':'large',
#          'ytick.labelsize':'large',
#          'figure.facecolor':'w'} 
# plt.rcParams.update(params)

lbg1 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/misc/LSST_Y4_CC_COSMOS_lbg_targets.fits'))
lbg2 = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/misc/LSST_Y4_RF_COSMOS_lbg_targets.fits'))
print(len(lbg1))
print(len(lbg2))

# mask = np.full(len(lbg1), True)
# plt.figure(figsize=(10, 10))
# plt.plot(lbg1['RA'][mask], lbg1['DEC'][mask], '.', ms=2, alpha=1)
# plt.grid(alpha=0.5)
# plt.gca().invert_xaxis()
# plt.show()

# mask = np.full(len(lbg2), True)
# plt.figure(figsize=(10, 10))
# plt.plot(lbg2['RA'][mask], lbg2['DEC'][mask], '.', ms=2, alpha=1)
# plt.grid(alpha=0.5)
# plt.gca().invert_xaxis()
# plt.show()

# https://github.com/rongpu/Python/blob/master/user_modules/match_coord.py
sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

# Remove duplicates
idx1, idx2, d2d, d_ra, d_dec = match_coord(lbg1['RA'], lbg1['DEC'], lbg2['RA'], lbg2['DEC'], search_radius=0.5, plot_q=True)

mask = np.full(len(lbg1), True)
mask[idx1] = False
lbg1 = lbg1[mask]
print(len(lbg1), np.sum(mask)/len(mask))

lbg1 = lbg1[['RA', 'DEC']]
lbg2 = lbg2[['RA', 'DEC']]
cat = vstack([lbg1, lbg2])
print(len(cat))

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/LSST_Y4_CC_RF_COSMOS_lbg_targets_unique.fits')
