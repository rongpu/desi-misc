# CCDs with ra, dec within 2.2 deg of the field centers
# Two components: deep field exposures and wide field exposures
# Deep field exposures must meet stringent quality cuts

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


deep_ra = np.array([150.1166])
deep_dec = np.array([2.2058])

###################### CCDs covering the deep fields ######################

ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', columns=['ra', 'dec']))
ccd_search_radius = 2.2*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, ccd['ra'], ccd['dec'], search_radius=ccd_search_radius, plot_q=False, keep_all_pairs=True)
idx = np.sort(idx2)
ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', rows=idx))
print(len(ccd))

# basic quality cuts
basic_quality = np.full(len(ccd), False)
mask = ccd['filter']!='Y'
basic_quality[mask] = (ccd['ccd_cuts'][mask]==0) | (ccd['ccd_cuts'][mask]==2**14)  # ignore DEPTH_CUT
mask = ccd['filter']=='Y'
ccd_cuts = ccd['ccd_cuts'][mask].copy()
ccd_cuts[ccd_cuts&2**1>0] -= 2**1  # ignore NOT_GRZ
ccd_cuts[ccd_cuts&2**14>0] -= 2**14  # ignore DEPTH_CUT
ccd_cuts[ccd_cuts&2**15>0] -= 2**15  # ignore TOO_MANY_BAD_CCDS
basic_quality[mask] = (ccd_cuts==0)

ccd = ccd[basic_quality]
print(len(ccd))

###################### Select deep exposures ######################

exp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures-medians.fits'))
print(len(exp))

columns = ['image_filename', 'camera', 'expnum', 'plver', 'procdate', 'plprocid', 'object', 'propid', 'filter', 'exptime', 'mjd_obs', 'airmass', 'ra_bore', 'dec_bore', 'zpt']
tmp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures.fits', columns=columns))
print(len(tmp))

exp = join(exp, tmp, keys='expnum')

mask = exp['n_good_ccds']>30
print(np.sum(mask)/len(mask))
exp = exp[mask]
print(len(exp))

fwhm = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-median-psfex-fwhm.fits'))
print(len(fwhm))
exp = join(exp, fwhm, join_type='inner', keys='expnum')
print(len(exp))

# Quality cuts for deep exposures
fwhm_limits = {'g':1.6, 'r':1.5, 'i':1.4, 'z':1.4, 'Y':2}
ccdskycounts_limits = {'g':2, 'r':5.5, 'i':15, 'z':30, 'Y':30}
zpt_limits = {'g':24.95, 'r':25.15, 'i':25.15, 'z':24.85, 'Y':23.7}

exp['good'] = False
for band in ['g', 'r', 'i', 'z', 'Y']:
    print(band)
    mask = exp['filter']==band
    print(np.sum(mask))
    exp['good'][mask] = (exp['median_psf_fwhm'][mask]*0.262<fwhm_limits[band])
    print(np.sum(exp['good'][mask]))
    exp['good'][mask] &= (exp['median_ccdskycounts'][mask]<ccdskycounts_limits[band])
    print(np.sum(exp['good'][mask]))
    exp['good'][mask] &= (exp['zpt'][mask]>zpt_limits[band])
    print(np.sum(exp['good'][mask]))
    # exp['good'][mask] &= (exp['exptime'][mask]>=60)
    # print(np.sum(exp['good'][mask]))

# sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
# from match_coord import match_coord
# exp_search_radius = 4.*3600.  # arcsec
# idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, exp['ra_bore'], exp['dec_bore'], search_radius=exp_search_radius, plot_q=False, keep_all_pairs=True)
# exp = exp[idx2]
# print('exp', len(exp))

###################### Select wide-field exposures ######################

des_wide = exp['propid']=='2012B-0001'
des_wide &= (np.char.startswith(exp['object'], 'DES survey hex')) | (np.char.startswith(exp['object'], 'DES wide hex'))
decals = exp['propid']=='2014B-0404'

mask = exp['good'].copy()
mask |= des_wide | decals
exp = exp[mask]
print('exp', len(exp))

###################### Get final CCD list ######################

mask = np.in1d(ccd['expnum'], exp['expnum'])
ccd = ccd[mask]
print('ccd', np.sum(mask))

exp = exp[['expnum', 'median_fwhm', 'median_ccdskycounts', 'median_psf_fwhm', 'n_ccds', 'n_good_ccds']]
ccd = join(ccd, exp, keys='expnum', join_type='left')
print('ccd', np.sum(mask))

ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-v4-cosmos-option-b.fits', overwrite=True)


