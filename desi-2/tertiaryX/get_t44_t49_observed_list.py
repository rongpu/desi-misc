# Get list of objects already observed in tertiary44 and tertiary49 so they will be excluded from tertiaryX

import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp


##########################################################################################################################################

fa = Table(fitsio.read('/global/cfs/cdirs/desi/survey/fiberassign/special/tertiary/0044/tertiary-targets-0044-assign.fits'))
print(len(fa))
fa = fa[fa['NASSIGN']>0]
print(len(fa))
print()

cat = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/matterhorn/zcatalog/v2/special/zpix-special-other.fits'))
tmp = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/matterhorn/zcatalog/v2/special/zpix-special-other-extra.fits'))
cat['TSNR2_LRG'] = tmp['TSNR2_LRG']
print(len(cat))

mask = np.in1d(cat['TARGETID'], fa['TARGETID'])
print(np.sum(mask))
cat = cat[mask]

# select science targets
mask = cat['OBJTYPE']=='TGT'
print('Science targets', np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]

mask = cat['COADD_FIBERSTATUS']==0
print('COADD_FIBERSTATUS', np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]

cat['EFFTIME_LRG'] = cat['TSNR2_LRG'] * 12.15
mask = cat['EFFTIME_LRG'] > 2100
print('EFFTIME_LRG', np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]

t44 = cat.copy()
t44 = t44[['TARGET_RA', 'TARGET_DEC', 'EFFTIME_LRG']]
t44['program'] = 't44'

##########################################################################################################################################

fa = Table(fitsio.read('/global/cfs/cdirs/desi/survey/fiberassign/special/tertiary/0049/tertiary-targets-0049-assign.fits'))
print(len(fa))
fa = fa[fa['NASSIGN']>0]
print(len(fa))
print()

cat = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/matterhorn/zcatalog/v2/special/zpix-special-other.fits'))
tmp = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/matterhorn/zcatalog/v2/special/zpix-special-other-extra.fits'))
cat['TSNR2_LRG'] = tmp['TSNR2_LRG']
print(len(cat))

mask = np.in1d(cat['TARGETID'], fa['TARGETID'])
print(np.sum(mask))
cat = cat[mask]

# select science targets
mask = cat['OBJTYPE']=='TGT'
print('Science targets', np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]

mask = cat['COADD_FIBERSTATUS']==0
print('COADD_FIBERSTATUS', np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]

cat['EFFTIME_LRG'] = cat['TSNR2_LRG'] * 12.15
mask = cat['EFFTIME_LRG'] > 2100
print('EFFTIME_LRG', np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]

t49 = cat.copy()
t49 = t49[['TARGET_RA', 'TARGET_DEC', 'EFFTIME_LRG']]
t49['program'] = 't49'

##########################################################################################################################################

obs = vstack([t44, t49])
print(len(obs))
obs.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/misc/observed_in_t44_t49.fits')
