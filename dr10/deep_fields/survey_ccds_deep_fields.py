# CCDs with ra, dec within 2.2 deg of the field centers

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

deep_ra = np.array([54.2743, 54.2743, 52.6484, 34.4757, 35.6645, 36.4500, 42.8200, 41.1944, 7.8744, 9.5000, 150.1166])
deep_dec = np.array([-27.1116, -29.0884, -28.1000, -4.9295, -6.4121, -4.6000, 0.0000, -0.9884, -43.0096, -43.9980, 2.2058])
field_names = ['C1', 'C2', 'C3', 'X1', 'X2', 'X3', 'S1', 'S2', 'E1', 'E2', 'COSMOS']

# bad exposures found in visual inspection
bad_expids = [938452, 964237, 964238, 964239, 980290, 980291, 986155, 1002854, 1006183, 1006184]

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

mask = np.in1d(exp['expnum'], bad_expids)
exp = exp[~mask]
print(len(exp))

psfex = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-psfex-fwhm.fits'))
_, idx = np.unique(psfex['expnum'], return_index=True)
psfex_exp = psfex[idx].copy()
exp = join(exp, psfex_exp, join_type='inner', keys='expnum')
print(len(exp))

# Quality cuts for deep exposures
seeing_limits = {'g':1.65, 'r':1.55, 'i':1.4, 'z':1.4, 'Y':2.0}
ccdskycounts_limits = {'g':2, 'r':5.5, 'i':15, 'z':30, 'Y':30}
zpt_limits = {'g':24.9, 'r':25.15, 'i':25.15, 'z':24.85, 'Y':23.7}

exp['quality_deep'] = False
exp['quality_wide'] = False
for band in ['g', 'r', 'i', 'z', 'Y']:
    print(band)
    mask = exp['filter']==band
    print(np.sum(mask))
    exp['quality_deep'][mask] = exp['median_psf_fwhm'][mask]*0.262<seeing_limits[band]
    exp['quality_wide'][mask] = exp['median_psf_fwhm'][mask]*0.262<seeing_limits[band]*1.2
    print(np.sum(exp['quality_deep'][mask]))
    exp['quality_deep'][mask] &= exp['median_ccdskycounts'][mask]<ccdskycounts_limits[band]
    exp['quality_wide'][mask] &= exp['median_ccdskycounts'][mask]<ccdskycounts_limits[band]*1.5
    print(np.sum(exp['quality_deep'][mask]))
    exp['quality_deep'][mask] &= exp['zpt'][mask]>zpt_limits[band]
    exp['quality_wide'][mask] &= exp['zpt'][mask]>zpt_limits[band]-0.2
    # exp['quality_deep'][mask] &= exp['exptime'][mask]>=60
    # exp['quality_wide'][mask] &= exp['exptime'][mask]>=60
    print(np.sum(exp['quality_deep'][mask]))

# Add i-band off-center exposures in COSMOS to fill in the gaps
exp_search_radius = 4.*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord([150.1166], [2.2058], exp['ra_bore'], exp['dec_bore'], search_radius=exp_search_radius, plot_q=False, keep_all_pairs=True)
exp['off_center'] = False
exp['off_center'][idx2[d2d>1.5*3600]] = True
exp['off_center'] &= exp['filter']=='i'

des_wide = exp['propid']=='2012B-0001'
des_wide &= (np.char.startswith(exp['object'], 'DES survey hex')) | (np.char.startswith(exp['object'], 'DES wide hex'))
decals = exp['propid']=='2014B-0404'

mask = exp['quality_deep'].copy()
mask |= (des_wide | decals | exp['off_center']) & exp['quality_wide']
exp = exp[mask]
print('exp', len(exp))

# Require fringe fitting in DR9 (for normalization)


###################### CCDs covering the deep fields ######################

ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', columns=['ra', 'dec']))
ccd_search_radius = 2.2*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, ccd['ra'], ccd['dec'], search_radius=ccd_search_radius, plot_q=False, keep_all_pairs=True)
idx = np.sort(idx2)
ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', rows=idx))
print('ccd', len(ccd))

mask = np.in1d(ccd['expnum'], exp['expnum'])
ccd = ccd[mask]
print('ccd', np.sum(mask))

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
print('ccd', len(ccd))

# Remove exposures with artifacts
mask = np.char.startswith(ccd['image_filename'], 'decam/CP/V4.8.2a/CP20131128')  # many of these exposures look bad
ccd = ccd[~mask]
print(len(ccd))

# require PLVER>=V4.8.1 except for Y band
mask = (ccd['plver']>='V4.8.1') | (ccd['filter']=='Y')
ccd = ccd[mask]
print(len(ccd))

# Add PSFEx FWHM values
ccd['ccd_id_str'] = np.char.add(np.array(ccd['expnum']).astype(str), ccd['ccdname'])
psfex['ccd_id_str'] = np.char.add(np.array(psfex['expnum']).astype(str), psfex['ccdname'])
psfex = psfex[['ccd_id_str', 'psf_fwhm', 'moffat_alpha', 'moffat_beta', 'failure']]
psfex.rename_column('failure', 'moffat_failure')
ccd = join(ccd, psfex, keys='ccd_id_str', join_type='left')
print('ccd', len(ccd))
ccd.remove_column('ccd_id_str')

# Remove CCDs with extreme FWHM (e.g., bad PSFEx models)
mask = (ccd['psf_fwhm']*0.262 > 0.6) & (ccd['psf_fwhm']*0.262 < 2.4)
ccd = ccd[mask]
print(len(ccd))

###################### Get final CCD list ######################

exp = exp[['expnum', 'median_fwhm', 'median_ccdskycounts', 'median_psf_fwhm', 'n_ccds', 'n_good_ccds', 'quality_deep', 'quality_wide']]
ccd = join(ccd, exp, keys='expnum', join_type='left')
print('ccd', len(ccd))

ccd['efftime'] = 10**(0.4*ccd['ccdzpt']-9) * ccd['exptime'] / (ccd['median_ccdskycounts'] * ccd['psf_fwhm']**2)
# ccd['exposure_efftime'] = 10**(0.4*ccd['zpt']-9) * ccd['exptime'] / (ccd['median_ccdskycounts'] * ccd['median_psf_fwhm']**2)

ccd.rename_column('ccd_cuts', 'ccd_cuts_dr10')
ccd['ccd_cuts'] = 0

ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1.fits', overwrite=True)
