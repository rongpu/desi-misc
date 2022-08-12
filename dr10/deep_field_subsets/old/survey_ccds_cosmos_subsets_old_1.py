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

###################################################################################################################

ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v6.fits', columns=['ra', 'dec']))
search_radius = 1.3*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, ccd['ra'], ccd['dec'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
idx = np.sort(idx2)
ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v6.fits', rows=idx))
print('ccd', len(ccd))

ccd.rename_column('ccd_cuts', 'ccd_cuts_dr10')
ccd['ccd_cuts'] = 0

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

mask = ~np.in1d(ccd['ccd_id_str'], dr9['ccd_id_str'])
ccd = ccd[mask]

################################## Get DR9 Legacy Survey CCD list as reference ##################################

exp_dr9 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-trim.fits'))

mask = exp_dr9['ccd_cuts']==0
exp_dr9 = exp_dr9[mask]

mask = (exp_dr9['propid']=='2014B-0404') | (exp_dr9['propid']=='2012B-0001')
exp_dr9 = exp_dr9[mask]
print('exp_dr9', len(exp_dr9))

mask = (exp_dr9['dec']>-20)
exp_dr9 = exp_dr9[mask]
print('exp_dr9', len(exp_dr9))

exp_columns = ['expnum', 'median_fwhm', 'median_ccdskycounts', 'median_psf_fwhm', 'n_ccds', 'n_good_ccds', 'exptime', 'efftime']
exp_dr9 = join(exp_dr9, exp[exp_columns], keys='expnum', join_type='left')
print('exp_dr9', len(exp_dr9))

###################################################################################################################

# Only keep exposures within 0.2 deg of field centers
search_radius = 0.2*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, exp['ra_bore'], exp['dec_bore'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
exp_cosmos = exp[idx2]
print('exp_cosmos', len(exp_cosmos))

# require PLVER>=V4.8.1
mask = (exp_cosmos['plver']>='V4.8.1')
exp_cosmos = exp_cosmos[mask]
print('exp_cosmos', len(exp_cosmos))

############################## quality cuts for the random subsets ##############################

mask = exp_cosmos['n_good_ccds']>30
print(np.sum(mask)/len(mask))
exp_cosmos = exp_cosmos[mask]
print(len(exp_cosmos))

mask = np.in1d(exp_cosmos['expnum'], bad_expids)
exp_cosmos = exp_cosmos[~mask]
print(len(exp_cosmos))

# Quality cuts for deep exposures
seeing_limits = {'g': 1.65, 'r': 1.55, 'i': 1.4, 'z': 1.4, 'Y': 2.0}
ccdskycounts_limits = {'g': 2, 'r': 5.5, 'i': 15, 'z': 30, 'Y': 30}
zpt_limits = {'g': 24.9, 'r': 25.15, 'i': 25.15, 'z': 24.85, 'Y': 23.7}

exptime_limits = {'g': 30, 'r': 30, 'z': 50}
efftime_limits = {'g': 7, 'r': 2, 'z': 0.5}

exp_cosmos['quality'] = False
for band in ['g', 'r', 'z']:
    print(band)
    mask = exp_cosmos['filter']==band
    print(np.sum(mask))
    exp_cosmos['quality'][mask] = exp_cosmos['median_psf_fwhm'][mask]*0.262<seeing_limits[band]*1.6
    exp_cosmos['quality'][mask] &= exp_cosmos['median_ccdskycounts'][mask]<ccdskycounts_limits[band]*1.5
    exp_cosmos['quality'][mask] &= exp_cosmos['zpt'][mask]>zpt_limits[band]-0.2
    exp_cosmos['quality'][mask] &= exp_cosmos['exptime'][mask]>=exptime_limits[band]
    exp_cosmos['quality'][mask] &= exp_cosmos['efftime'][mask]>=efftime_limits[band]

mask = exp_cosmos['quality'].copy()
exp_cosmos = exp_cosmos[mask]
print('exp_cosmos', len(exp_cosmos))

# Downselect the exposures
np.random.seed(9371)
expnum_list = []
for band in ['g', 'r', 'z']:
    idx = np.where(exp_dr9['filter']==band)[0]
    idx = np.random.choice(idx, size=50, replace=False)
    for index in idx:
        mask = exp_cosmos['filter']==band
        mask &= ~np.in1d(exp_cosmos['expnum'], expnum_list)
        index_best = np.argmin(np.abs(exp_cosmos['efftime'][mask]-exp_dr9['efftime'][index]))
        expnum_list.append(exp_cosmos['expnum'][mask][index_best])

mask = np.in1d(exp_cosmos['expnum'], expnum_list)
exp_cosmos = exp_cosmos[mask]
print('exp_cosmos', len(exp_cosmos))

########## CCDs list ##########

mask = np.in1d(ccd['expnum'], exp_cosmos['expnum'])
ccd = ccd[mask]
print('ccd', np.sum(mask))

# basic quality cuts
ccd_cuts = ccd['ccd_cuts_dr10'].copy()
ccd_cuts[ccd_cuts&2**13>0] -= 2**13  # ignore EARLY_DECAM
ccd_cuts[ccd_cuts&2**14>0] -= 2**14  # ignore DEPTH_CUT
mask = ccd_cuts==0
ccd = ccd[mask]
print('ccd', len(ccd))

# Remove exposures with artifacts
mask = np.char.startswith(ccd['image_filename'], 'decam/CP/V4.8.2a/CP20131128')  # many of these exposures look bad
ccd = ccd[~mask]
print('ccd', len(ccd))

# require PLVER>=V4.8.1
mask = (ccd['plver']>='V4.8.1')
ccd = ccd[mask]
print('ccd', len(ccd))

# Remove CCDs with extreme FWHM (e.g., bad PSFEx models)
mask = (ccd['psf_fwhm']*0.262 > 0.6) & (ccd['psf_fwhm']*0.262 < 2.4)
ccd = ccd[mask]
print('ccd', len(ccd))

########## Get final CCD list ##########

################################## random subsets ##################################

seeds = {'g': 188, 'r': 39, 'z': 88}

idx_split = {}
exp1 = {}
min_subset = 100000
for band in ['g', 'r', 'z']:

    np.random.seed(seeds[band])

    # Shuffle
    _, idx = np.unique(ccd['expnum'], return_index=True)
    exp1[band] = ccd[idx].copy()
    idx = np.random.choice(len(exp1[band]), size=len(exp1[band]), replace=False)
    exp1[band] = exp1[band][idx]

    idx = np.where(exp1[band]['filter']==band)[0]
    if band=='r':
        split_idx = np.cumsum(np.random.choice([1, 2, 3], p=[0.4, 0.5, 0.1], size=1000))
    else:
        split_idx = np.cumsum(np.random.choice([2, 3, 4, 5, 6], p=[0.1, 0.2, 0.25, 0.35, 0.1], size=1000))
    split_idx = split_idx[split_idx<len(idx)]
    if (split_idx[-1]!=len(idx)-1) and (split_idx[-1]<len(idx)-3):
        split_idx = np.append(split_idx, len(idx)-1)
    split = np.array_split(idx, split_idx)
    efftime_list = np.zeros(len(split))
    for index in range(len(split)):
        mask = np.in1d(ccd['expnum'], exp1[band]['expnum'][split[index]])
        mask &= ccd['inner']
        efftime_list[index] = np.sum(ccd['efftime'][mask])
    sort_subsets = np.argsort(efftime_list)
    idx_split[band] = [split[index] for index in sort_subsets]
    min_subset = np.minimum(min_subset, len(idx_split[band]))
    print('split', band, len(idx_split[band]))

np.random.seed(34399)

print('min_subset', min_subset)
for band in ['g', 'r', 'z']:
    print(band, len(idx_split[band]))
    while len(idx_split[band])>min_subset:
        if band=='r':
            efftime_list = np.zeros(len(idx_split[band]))
            for index in range(len(idx_split[band])):
                mask = np.in1d(ccd['expnum'], exp1[band]['expnum'][idx_split[band][index]])
                mask &= ccd['inner']
                efftime_list[index] = np.sum(ccd['efftime'][mask])
            if efftime_list[0]<1100:
                index_pop = 0
            elif efftime_list[-1]>6000:
                index_pop = -1
            else:
                index_pop = np.random.choice(len(idx_split[band]))
        elif band=='z':
            index_pop = 0
        else:
            index_pop = np.random.choice(len(idx_split[band]))
        idx_split[band].pop(index_pop)
    print(len(idx_split[band]))

for index in range(min_subset):
    mask = np.full(len(ccd), False)
    for band in ['g', 'r', 'z']:
        mask |= np.in1d(ccd['expnum'], exp1[band]['expnum'][idx_split[band][index]])
    ccd_subset = ccd[mask].copy()
    # ccd_subset.remove_columns('inner')
    ccd_subset.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/survey-ccds-dr10-v6-subset-{}.fits'.format(index), overwrite=True)

# ccd_dr9.remove_column('ccd_id_str')
# ccd_dr9.remove_columns('inner')
ccd_dr9.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/survey-ccds-dr10-v6-subset-dr9-ccds.fits', overwrite=True)
