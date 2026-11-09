from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

params = {'legend.fontsize': 'medium',
          'axes.labelsize': 'medium',
          'axes.titlesize': 'medium',
          'xtick.labelsize': 'medium',
          'ytick.labelsize': 'medium',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

fns = [
    '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary44/tertiary44_lae_truth-20261107.fits',
    '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary47/catalogs/tertiary47_ibis_lae_truth-20251104.fits',
    '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/catalogs/tertiary49_ibis_lae_truth-20260127-add_target_info-add_lyaflux-newcol.fits',
    '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/catalogs/tertiary54_ibis_filler_lae_truth-20260413-newcol.fits',
    '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/catalogs/tertiary54_ibis_primary_lae_truth-20260413-newcol.fits',
]

cat_stack = []
for fn in fns:
    print(fn)
    cat = Table(fitsio.read(fn))
    print(len(cat))
    
    if 'tertiary44_lae_truth-20261107.fits' in fn:
        mask = (cat['TERTIARY_TARGET']=='HIZ_IBIS') | (cat['TERTIARY_TARGET']=='HIZ1H_IBIS')
        cat = cat[mask]
        print('tertiary44', len(cat))
    
    mask = cat['COADD_FIBERSTATUS']==0
    print(np.sum(mask), np.sum(mask)/len(mask))
    cat = cat[mask]

    cat['EFFTIME_LRG'] = cat['TSNR2_LRG'] * 12.15
    mask = cat['EFFTIME_LRG'] > 1800 * 0.8
    print(np.sum(mask), np.sum(mask)/len(mask))
    cat = cat[mask]
    print(len(cat))
    cat.remove_column('TSNR2_LRG')
    
    mask = cat['lae'].copy()
    print(cat['Z'][mask].min(), cat['Z'][mask].max())
    
    if np.sum(cat['TARGET_RA']<100)>0 and np.sum(cat['TARGET_RA']>100)==0:
        cat['field'] = 'XMM-LSS'
    elif np.sum(cat['TARGET_RA']>100)>0 and np.sum(cat['TARGET_RA']<100)==0:
        cat['field'] = 'COSMOS'
    else:
        raise ValueError

    cat['program'] = 'tertiary'+fn[fn.find('tertiary')+8:fn.find('tertiary')+10]
    cat = cat[['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'EFFTIME_LRG', 'program', 'field', 'Z', 'lae']]
    cat_stack.append(cat)

cat = vstack(cat_stack)
print(len(cat))

################################### Remove duplicates ###################################
sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord, match_self

score = cat['EFFTIME_LRG'].copy()
n_duplicates = 1

while n_duplicates>0:
    n_duplicates, idx1, idx2 = match_self(cat['TARGET_RA'], cat['TARGET_DEC'], search_radius=1., return_indices=True, plot_q=False)
    if len(idx1)==0:
        continue
    mask_remove1 = (score[idx1] < score[idx2]) | (score[idx1]==0)
    mask_remove2 = (score[idx2] < score[idx1]) | (score[idx2]==0)
    idx_remove = np.unique(np.concatenate([idx1[mask_remove1], idx2[mask_remove2]]))
    print(len(idx_remove), len(idx_remove)/len(idx1))
    mask = np.full(len(cat), True)
    mask[idx_remove] = False
    print(np.sum(mask), np.sum(mask)/len(mask))
    cat = cat[mask]
    score = score[mask]

print()
print(len(cat))

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/misc/lae_truth/lae_truth_combined_20261109.fits')

# plt.plot(cat['TARGET_RA'], cat['TARGET_DEC'], '.', ms=1)

tmp = Table()
tmp['program'], tmp['count'] = np.unique(cat['program'], return_counts=True)
tmp.pprint_all()

progs, idx, inv, counts = np.unique(
    cat['program'], return_index=True, return_inverse=True, return_counts=True
)

tmp = Table()
tmp['program'] = progs
tmp['field'] = cat['field'][idx]
tmp['n_obs'] = counts
tmp['n_lae'] = np.bincount(inv, weights=np.asarray(cat['lae'], dtype=float)).astype(int)
tmp['lae_frac (%)'] = tmp['n_lae']/tmp['n_obs']*100
tmp['lae_frac (%)'].format = '%.1f'
tmp.pprint_all()
print()

print('Total n_obs:', tmp['n_obs'].sum())
print('Total n_lae:', tmp['n_lae'].sum())

