# Create a combined survey-ccd list
# Only keep the deep and wide exposures (exclude Pyxis, Pal14, NGC7492, etc.)
# Only keep the good (e.g., ccd_cuts==0) exposures

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

ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']

############## COSMOS ##############

dirname = '/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/'
fns = [
'survey-ccds-ibis-3.fits',
'survey-ccds-ibis-4.fits',
'survey-ccds-ibis-5.fits',
'survey-ccds-ibis-6.fits',
'survey-ccds-ibis-7.fits',
'survey-ccds-ibis-8.fits',
'survey-ccds-ibis-9.fits',
]

ccd_stack = []

for fn in fns:

    ccd = Table(fitsio.read(os.path.join(dirname, fn)))
    print(len(ccd))

    ccd_stack.append(ccd)

ccd = vstack(ccd_stack)
print(len(ccd), len(np.unique(ccd['expnum'])))

mask = np.in1d(ccd['filter'], ibis_filters)
ccd = ccd[mask]
print(len(ccd), np.sum(mask)/len(mask))

mask = ccd['ccd_cuts']==0
ccd = ccd[mask]
print(len(ccd), np.sum(mask)/len(mask))

mask = ccd['exptime']>100
ccd = ccd[mask]
print(len(ccd), np.sum(mask)/len(mask))

print('\nNon-standard CCDs')
tmp = Table()
tmp['object'], tmp['count'] = np.unique(ccd['object'], return_counts=True)
nonstandard = np.full(len(ccd), False)
for index in range(len(tmp)):
    if (not tmp['object'][index].startswith('IBIS_deep')) and (not tmp['object'][index].startswith('IBIS_wide')):
        # print(tmp['object'][index], tmp['count'][index])
        mask = ccd['object']==tmp['object'][index]
        nonstandard[mask]=True
print('\nNumber of non-standard CCDs: {} ({:.2f}%)'.format(np.sum(nonstandard), 100*np.sum(nonstandard/len(nonstandard))))
ccd = ccd[~nonstandard]
print(len(ccd), np.sum(~nonstandard)/len(nonstandard))

tmp = Table()
tmp['object'], tmp['count'] = np.unique(ccd['object'], return_counts=True)
deep = np.full(len(ccd), False)
wide = np.full(len(ccd), False)
for index in range(len(tmp)):
    mask = ccd['object']==tmp['object'][index]
    if tmp['object'][index].startswith('IBIS_deep'):
        deep[mask]=True
    elif tmp['object'][index].startswith('IBIS_wide'):
        wide[mask]=True
    else:
        raise ValueError
ccd['ibis_deep'] = deep.copy()
ccd['ibis_wide'] = wide.copy()

ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/misc/survey-ccds-ibis-deep-and-wide-3-9.fits')

