from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


cat = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/catalogs/cosmos.fits'))
print(len(cat))
cat['original_index'] = np.arange(len(cat))

hsc = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-dud-no_mag_limit-reduced.fits', columns=['ra']))
print(len(hsc))
idx = np.where((hsc['ra']>140) & (hsc['ra']<160))[0]
hsc = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-dud-no_mag_limit-reduced.fits', rows=idx))
print(len(hsc))

idx1, idx2, d2d, d_ra, d_dec = match_coord(hsc['ra'], hsc['dec'], cat['ra'], cat['dec'], search_radius=1., plot_q=True)
hsc = hsc[idx1]
cat = cat[idx2]

cat.write('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/catalogs/hsc_pdr3_matched/tractor_matched.fits')
hsc.write('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10-deep/cosmos/catalogs/hsc_pdr3_matched/hsc_matched.fits')
