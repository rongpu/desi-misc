from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio


surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-defringed.fits'
ccd_all = Table(fitsio.read(surveyccd_path))
print('ccd all', len(ccd_all))


fns = glob.glob('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/tmp/*.fits')

ccd = []
for fn in fns:
    tmp = Table(fitsio.read(fn))
    ccd.append(tmp)
ccd = vstack(ccd)

ccd_all['ccd_id_str'] = np.char.add(np.array(ccd_all['expnum']).astype(str), ccd_all['ccdname'])
ccd['ccd_id_str'] = np.char.add(np.array(ccd['expnum']).astype(str), ccd['ccdname'])

print(len(ccd), len(ccd['ccd_id_str']))

_, idx = np.unique(ccd['ccd_id_str'], return_index=True)
if len(idx)>len(ccd):
    print('duplicates', len(idx), len(ccd))
    ccd = ccd[idx]

mask = np.in1d(ccd_all['ccd_id_str'], ccd['ccd_id_str'])
ccd_all = ccd_all[mask]
if not len(ccd_all)==len(ccd):
    raise ValueError

# Restore order
if len(ccd_all)!=len(ccd) or not np.all(np.unique(ccd_all['ccd_id_str'])==np.unique(ccd['ccd_id_str'])):
    raise ValueError('ccd_all and ccd have different ccd_id_str list')
ccd_all_reverse_sort = np.array(ccd_all['ccd_id_str']).argsort().argsort()
ccd = ccd[np.argsort(ccd['ccd_id_str'])[ccd_all_reverse_sort]]
if not np.all(ccd['ccd_id_str']==ccd_all['ccd_id_str']):
    raise ValueError

print(len(ccd), len(ccd['ccd_id_str']))

ccd.remove_column('ccd_id_str')

ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-defringed-final.fits', overwrite=True)
