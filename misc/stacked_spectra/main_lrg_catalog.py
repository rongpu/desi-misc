from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/spectro/everest/main_cumulative_lrg.fits'))
cat1 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/spectro/everest/sv3_cumulative_lrg.fits'))
cat1 = cat1[cat1['main_lrg']]
cat1.remove_columns(['main_lrg', 'SV3_DESI_TARGET', 'SV3_BGS_TARGET'])
cat = vstack([cat, cat1])

# cat = Table(fitsio.read('/Users/rongpu/Documents/Data/desi_data/everest/main_cumulative_lrg.fits'))
cat['EFFTIME_ELG'] = 8.60 * cat['TSNR2_ELG']
cat['EFFTIME_LRG'] = 12.33 * cat['TSNR2_LRG']

# Remove FIBERSTATUS!=0 fibers
mask = cat['COADD_FIBERSTATUS']==0
print('FIBERSTATUS',np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Remove "no data" fibers
mask = cat['ZWARN'] & 2**9==0
print('No data', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Apply LRG mask
mask = cat['lrg_mask']==0
print('LRG mask', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Remove QSO targets
mask = cat['DESI_TARGET'] & 2**2 ==0
print('Remove QSO targets', np.sum(mask), np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

# Require a minimum depth
min_depth = 800.
mask = cat['EFFTIME_LRG']>min_depth
print('Min depth', np.sum(mask), np.sum(~mask), np.sum(mask)/len(mask))
cat = cat[mask]

# Julien's bad fibers list
bad_fibers = np.array(Table.read('/global/cfs/cdirs/desi/users/rongpu/spectro/everest/misc/badfibers.csv')['FIBER'])
# bad_fibers = np.array(Table.read('/Users/rongpu/Documents/Data/desi_data/everest/misc/badfibers.csv')['FIBER'])
bad_fibers = np.append(bad_fibers, np.arange(2663, 2674+1))  # fibers affected by the CCD z5 defect
bad_fibers = np.append(bad_fibers, [3402, 3429])  # "swapped" fibers
bad_fibers = np.unique(bad_fibers)
print(len(bad_fibers), 'bad fibers')
mask_bad = np.in1d(cat['FIBER'], bad_fibers)
print('Bad fibers', np.sum(~mask_bad), np.sum(mask_bad), np.sum(mask_bad)/len(mask_bad))
cat = cat[~mask_bad]

# Remove duplidates keeping the higher EFFTIME objects
print(len(cat), len(np.unique(cat['TARGETID'])), len(cat)-len(np.unique(cat['TARGETID'])))
cat.sort('EFFTIME_LRG', reverse=True)
_, idx_keep = np.unique(cat['TARGETID'], return_index=True)
cat = cat[idx_keep]
print(len(cat), len(np.unique(cat['TARGETID'])), len(cat)-len(np.unique(cat['TARGETID'])))

print(len(cat))

# Custom DELTACHI2 vs z cut
d = (10**(3 - 3.5*cat['Z']))
mask_remove = (d>30) & (cat['DELTACHI2']<30)
mask_remove |= (d<30) & (cat['DELTACHI2']<d)
mask_remove |= (cat['DELTACHI2']<10)
mask_quality = cat['ZWARN']==0
mask_quality &= cat['Z']<1.4
mask_quality &= (~mask_remove)

print(np.sum(~mask_quality)/len(mask_quality))
cat = cat[mask_quality]
print(len(cat))

# Remove stars and "QSO"s
mask_star = (cat['SPECTYPE']=='STAR') | (cat['Z']<0.0003)
print('Star', np.sum(mask_star), np.sum(mask_star)/len(mask_star))
cat = cat[~mask_star]
print(len(cat))
mask_qso = cat['SPECTYPE']=='QSO'
print('Star', np.sum(mask_qso), np.sum(mask_qso)/len(mask_qso))
cat = cat[~mask_qso]
print(len(cat))

# redshift cut
mask = (cat['Z']>0.4) & (cat['Z']<1.0)
print('Redshift', np.sum(~mask), np.sum(~mask)/len(mask))
cat = cat[mask]

cat['thrunight'] = '00000000'
for tileid in np.unique(cat['TILEID']):
    dirname = glob.glob('/global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/{}/*'.format(str(tileid)))
    if len(dirname)!=1:
        raise ValueError
    dirname = dirname[0]
    thrunight = os.path.basename(dirname)
    mask = cat['TILEID']==tileid
    cat['thrunight'][mask] = thrunight

cat.write('/global/cfs/cdirs/desi/users/rongpu/tmp/lrgs_for_stacking.fits', overwrite=True)
