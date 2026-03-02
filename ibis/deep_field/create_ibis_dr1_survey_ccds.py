# Create survey-ccd list for IBIS DR1 production

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

############## XMM-LSS ##############

fn = '/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/survey-ccds-ibis-3.fits'

ccd = Table(fitsio.read(fn))
print(len(ccd))

mask = np.in1d(ccd['filter'], ibis_filters)
ccd = ccd[mask]
print(len(ccd))

mask = ccd['ccd_cuts']==0
ccd = ccd[mask]
print(len(ccd))

mask = ccd['exptime']>100
ccd = ccd[mask]
print(len(ccd))

ra, dec = 35.7, -4.8
ramin, ramax, decmin, decmax = ra-2.5, ra+2.5, dec-2.5, dec+2.5
mask = (ccd['ra_bore']>ramin) & (ccd['ra_bore']<ramax) & (ccd['dec_bore']>decmin) & (ccd['dec_bore']<decmax)
# plt.plot(ccd['ra_bore'][mask], ccd['dec_bore'][mask], '.', alpha=0.5)
# plt.axis([ramax, ramin, decmin, decmax])
# plt.title('XMM {} CCDs'.format(np.sum(mask)))
# plt.show()

ccd = ccd[mask]
print(len(ccd))

tmp = Table()
tmp['object'], tmp['count'] = np.unique(ccd['object'], return_counts=True)
nonstandard = False
for index in range(len(tmp)):
    if not tmp['object'][index].startswith('IBIS_deep'):
        print(tmp['object'][index], tmp['count'][index])
        mask = ccd['object']==tmp['object'][index]
        nonstandard[mask]=True
print('non-standard CCDs:', np.sum(nonstandard))

xmm = ccd.copy()
print('XMM', len(xmm))

############## COSMOS ##############

fns = ['/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/survey-ccds-ibis-3.fits',
       '/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/survey-ccds-ibis-4.fits',
       '/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/survey-ccds-ibis-5.fits']

ccd_stack = []

for fn in fns:

    ccd = Table(fitsio.read(fn))
    print(len(ccd))

    mask = np.in1d(ccd['filter'], ibis_filters)
    ccd = ccd[mask]
    print(len(ccd))

    mask = ccd['ccd_cuts']==0
    ccd = ccd[mask]
    print(len(ccd))

    mask = ccd['exptime']>100
    ccd = ccd[mask]
    print(len(ccd))

    ccd_stack.append(ccd)

ccd = vstack(ccd_stack)
    
ra, dec = 150.1, 2.2
ramin, ramax, decmin, decmax = ra-2.5, ra+2.5, dec-2.5, dec+2.5
mask = (ccd['ra_bore']>ramin) & (ccd['ra_bore']<ramax) & (ccd['dec_bore']>decmin) & (ccd['dec_bore']<decmax)
# plt.plot(ccd['ra_bore'][mask], ccd['dec_bore'][mask], '.', alpha=0.5)
# plt.axis([ramax, ramin, decmin, decmax])
# plt.title('COSMOS {} CCDs'.format(np.sum(mask)))
# plt.show()

ccd = ccd[mask]
print(len(ccd))

tmp = Table()
tmp['object'], tmp['count'] = np.unique(ccd['object'], return_counts=True)
nonstandard = False
for index in range(len(tmp)):
    if not tmp['object'][index].startswith('IBIS_deep'):
        print(tmp['object'][index], tmp['count'][index])
        mask = ccd['object']==tmp['object'][index]
        nonstandard[mask]=True
print('non-standard CCDs:', np.sum(nonstandard))

cosmos = ccd.copy()
print('COSMOS', len(cosmos))

############## Combined ##############

ccd = vstack([xmm, cosmos])
print('Combined', len(ccd))

ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/survey-ccds-ibis-dr1.fits')

