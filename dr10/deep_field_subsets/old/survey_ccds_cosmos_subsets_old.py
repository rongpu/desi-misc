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

# bad exposures found in visual inspection
bad_expids = [938452, 964237, 964238, 964239, 980290, 980291, 986155, 1002854, 1006183, 1006184]

exp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures-medians.fits'))
print(len(exp))

columns = ['image_filename', 'camera', 'expnum', 'plver', 'procdate', 'plprocid', 'object', 'propid', 'filter', 'exptime', 'mjd_obs', 'airmass', 'ra_bore', 'dec_bore', 'zpt']
tmp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures.fits', columns=columns))
print(len(tmp))
exp = join(exp, tmp, keys='expnum')

psfex = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-psfex-fwhm.fits'))
_, idx = np.unique(psfex['expnum'], return_index=True)
psfex_exp = psfex[idx].copy()
exp = join(exp, psfex_exp, join_type='inner', keys='expnum')
print(len(exp))

exp['efftime'] = 10**(0.4*exp['zpt']-9) * exp['exptime'] / (exp['median_ccdskycounts'] * exp['median_psf_fwhm']**2)

ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v6.fits', columns=['ra', 'dec']))
search_radius = 1.3*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, ccd['ra'], ccd['dec'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
idx = np.sort(idx2)
ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v6.fits', rows=idx))
print('ccd', len(ccd))

exp_columns = ['expnum', 'median_fwhm', 'median_ccdskycounts', 'median_psf_fwhm', 'n_ccds', 'n_good_ccds']
ccd = join(ccd, exp[exp_columns], keys='expnum', join_type='left')
print('ccd', len(ccd))

# Add PSFEx FWHM values
ccd['ccd_id_str'] = np.char.add(np.array(ccd['expnum']).astype(str), ccd['ccdname'])
psfex['ccd_id_str'] = np.char.add(np.array(psfex['expnum']).astype(str), psfex['ccdname'])
psfex = psfex[['ccd_id_str', 'psf_fwhm', 'moffat_alpha', 'moffat_beta', 'failure']]
psfex.rename_column('failure', 'moffat_failure')
ccd = join(ccd, psfex, keys='ccd_id_str', join_type='left')
print('ccd', len(ccd))

# Flag the inner 1 deg CCDs
search_radius = 1.0*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, ccd['ra'], ccd['dec'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
ccd['inner'] = False
ccd['inner'][idx2] = True

ccd['ccd_id_str'] = np.char.add(np.array(ccd['expnum']).astype(str), ccd['ccdname'])

ccd['efftime'] = 10**(0.4*ccd['ccdzpt']-9) * ccd['exptime'] / (ccd['median_ccdskycounts'] * ccd['psf_fwhm']**2)
# ccd['exposure_efftime'] = 10**(0.4*ccd['zpt']-9) * ccd['exptime'] / (ccd['median_ccdskycounts'] * ccd['median_psf_fwhm']**2)

################################## subset with identical CCDs as in DR9 ##################################

dr9 = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz', columns=['expnum', 'ccdname', 'ra', 'dec', 'ccd_cuts']))

mask = dr9['ccd_cuts']==0
dr9 = dr9[mask]

dr9['ccd_id_str'] = np.char.add(np.array(dr9['expnum']).astype(str), dr9['ccdname'])

search_radius = 1.3*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, dr9['ra'], dr9['dec'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
idx = np.sort(idx2)
dr9 = dr9[idx]
print(len(dr9))

mask = np.in1d(dr9['ccd_id_str'], ccd['ccd_id_str'])
if np.sum(mask)!=len(dr9):
    raise ValueError

mask = np.in1d(ccd['ccd_id_str'], dr9['ccd_id_str'])
ccd_dr9 = ccd[mask].copy()
ccd_dr9.remove_column('ccd_id_str')
# ccd_dr9.remove_columns('inner')

mask = ~np.in1d(ccd['ccd_id_str'], dr9['ccd_id_str'])
ccd = ccd[mask]

################################## quality cuts for the random subsets ##################################

########## Select exposures ##########

mask = exp['n_good_ccds']>30
print(np.sum(mask)/len(mask))
exp = exp[mask]
print(len(exp))

mask = np.in1d(exp['expnum'], bad_expids)
exp = exp[~mask]
print(len(exp))

# Quality cuts for deep exposures
seeing_limits = {'g': 1.65, 'r': 1.55, 'i': 1.4, 'z': 1.4, 'Y': 2.0}
ccdskycounts_limits = {'g': 2, 'r': 5.5, 'i': 15, 'z': 30, 'Y': 30}
zpt_limits = {'g': 24.9, 'r': 25.15, 'i': 25.15, 'z': 24.85, 'Y': 23.7}

exptime_limits = {'g':30, 'r':30, 'z':50}
efftime_limits = {'g':7, 'r':2, 'z':0.5}

exp['quality'] = False
for band in ['g', 'r', 'z']:
    print(band)
    mask = exp['filter']==band
    print(np.sum(mask))
    exp['quality'][mask] = exp['median_psf_fwhm'][mask]*0.262<seeing_limits[band]*1.6
    exp['quality'][mask] &= exp['median_ccdskycounts'][mask]<ccdskycounts_limits[band]*1.5
    exp['quality'][mask] &= exp['zpt'][mask]>zpt_limits[band]-0.2
    exp['quality'][mask] &= exp['exptime'][mask]>=exptime_limits[band]
    exp['quality'][mask] &= exp['efftime'][mask]>=efftime_limits[band]

mask = exp['quality'].copy()
exp = exp[mask]
print('exp', len(exp))

########## CCDs covering the deep fields ##########

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

# # require PLVER>=V4.8.1 except for Y band
# mask = (ccd['plver']>='V4.8.1') | (ccd['filter']=='Y')
# ccd = ccd[mask]
# print(len(ccd))

# Remove CCDs with extreme FWHM (e.g., bad PSFEx models)
mask = (ccd['psf_fwhm']*0.262 > 0.6) & (ccd['psf_fwhm']*0.262 < 2.4)
ccd = ccd[mask]
print(len(ccd))

########## Get final CCD list ##########

ccd.rename_column('ccd_cuts', 'ccd_cuts_dr10')
ccd['ccd_cuts'] = 0

################################## random subsets ##################################

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx].copy()

# exp_all = exp.copy()

# Only keep exposures within 0.2 deg of field centers
search_radius = 0.2*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, exp['ra_bore'], exp['dec_bore'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
exp = exp[idx2]
print('exp', len(exp))

np.random.seed(82361)

# Downselect exposures
mask_keep = np.full(len(exp), True)
idx = np.where(exp['filter']=='r')[0]
idx_remove = np.random.choice(idx, size=len(idx)//2, replace=False)
mask_keep[idx_remove] = False
idx = np.where(exp['filter']=='z')[0]
idx_remove = np.random.choice(idx, size=int(len(idx)*0.7), replace=False)
mask_keep[idx_remove] = False
exp = exp[mask_keep]
print('exp', len(exp))

# Downselect exposures
mask_keep = np.full(len(exp), True)
idx = np.where((exp['filter']=='g') & (exp['efftime']>90))[0]
idx_remove = np.random.choice(idx, size=int(len(idx)*0.4), replace=False)
mask_keep[idx_remove] = False
idx = np.where((exp['filter']=='r') & (exp['efftime']>50))[0]
idx_remove = np.random.choice(idx, size=int(len(idx)*0.4), replace=False)
mask_keep[idx_remove] = False
exp = exp[mask_keep]
print('exp', len(exp))


# Shuffle
idx = np.random.choice(len(exp), size=len(exp), replace=False)
exp = exp[idx]

n_subset = 24

idx_split = {}
for band in ['g', 'r', 'z']:
    idx = np.where(exp['filter']==band)[0]
    split = np.array_split(idx, n_subset)
    efftime_list = np.zeros(len(split))
    for index in range(len(split)):
        mask = np.in1d(ccd['expnum'], exp['expnum'][split[index]])
        mask &= ccd['inner']
        efftime_list[index] = np.sum(ccd['efftime'][mask])
    sort_subsets = np.argsort(efftime_list)
    idx_split[band] = [split[index] for index in sort_subsets]

for index in range(n_subset):
    mask = np.full(len(ccd), False)
    for band in ['g', 'r', 'z']:
        mask |= np.in1d(ccd['expnum'], exp['expnum'][idx_split[band][index]])
    ccd_subset = ccd[mask].copy()
    # ccd_subset.remove_columns('inner')
    ccd_subset.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/survey-ccds-dr10-v6-subset-{}.fits'.format(index), overwrite=True)

ccd_dr9.remove_column('ccd_id_str')
# ccd_dr9.remove_columns('inner')
ccd_dr9.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/survey-ccds-dr10-v6-subset-dr9-ccds.fits', overwrite=True)
