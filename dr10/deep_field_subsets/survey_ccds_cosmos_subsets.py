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

exp_all = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures-medians.fits'))
print(len(exp_all))

columns = ['image_filename', 'camera', 'expnum', 'plver', 'procdate', 'plprocid', 'object', 'propid', 'filter', 'exptime', 'mjd_obs', 'airmass', 'ra_bore', 'dec_bore', 'zpt']
tmp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures.fits', columns=columns))
print(len(tmp))
exp_all = join(exp_all, tmp, keys='expnum')

psfex = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-psfex-fwhm.fits'))
_, idx = np.unique(psfex['expnum'], return_index=True)
psfex_exp = psfex[idx].copy()
exp_all = join(exp_all, psfex_exp, join_type='inner', keys='expnum')
print(len(exp_all))

exp_all['efftime'] = 10**(0.4*exp_all['zpt']-9) * exp_all['exptime'] / (exp_all['median_ccdskycounts'] * exp_all['median_psf_fwhm']**2)

################################## subset with identical CCDs as in DR9 ##################################

dr9 = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz', columns=['expnum', 'ccdname', 'ra', 'dec', 'ccd_cuts']))

search_radius = 2.*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, dr9['ra'], dr9['dec'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
idx = np.sort(idx2)
print(len(idx))

mask = dr9['ccd_cuts'][idx]==0
idx = idx[mask]

dr9 = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz', rows=idx))
dr9.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/cosmos/survey-ccds-decam-dr9-cosmos.fits', overwrite=True)

dr9['ccd_id_str'] = np.char.add(np.array(dr9['expnum']).astype(str), dr9['ccdname'])

ccd_ref = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v6.fits', columns=['expnum', 'ccdname']))
ccd_ref['ccd_id_str'] = np.char.add(np.array(ccd_ref['expnum']).astype(str), ccd_ref['ccdname'])

mask = np.in1d(dr9['ccd_id_str'], ccd_ref['ccd_id_str'])
if np.sum(mask)!=len(dr9):
    raise ValueError

mask = np.in1d(ccd_ref['ccd_id_str'], dr9['ccd_id_str'])
idx = np.where(mask)[0]
ccd_ref = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v6.fits', rows=idx))

# Flag the inner 1 deg CCDs
search_radius = 1.0*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, ccd_ref['ra'], ccd_ref['dec'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
ccd_ref['inner'] = False
ccd_ref['inner'][idx2] = True

exp_columns = ['expnum', 'median_fwhm', 'median_ccdskycounts', 'median_psf_fwhm', 'efftime']
ccd_ref = join(ccd_ref, exp_all[exp_columns], keys='expnum', join_type='left')

# ccd_ref.remove_column('ccd_id_str')
# ccd_ref.remove_columns('inner')
ccd_ref.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/cosmos/survey-ccds-dr10-v6-dr9-ccds.fits', overwrite=True)

################################## Deep field CCD list ##################################

ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-defringed.fits'))
ccd['ccd_id_str'] = np.char.add(np.array(ccd['expnum']).astype(str), ccd['ccdname'])

_, idx = np.unique(ccd['expnum'], return_index=True)
exp_cosmos = ccd[idx].copy()

mask = ~np.in1d(exp_cosmos['expnum'], ccd_ref['expnum'])
exp_cosmos = exp_cosmos[mask]

# Only keep exposures within 0.2 deg of field centers
search_radius = 0.2*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, exp_cosmos['ra_bore'], exp_cosmos['dec_bore'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
exp_cosmos = exp_cosmos[idx2]
print('exp_cosmos', len(exp_cosmos))

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

exp_columns = ['expnum', 'median_fwhm', 'median_ccdskycounts', 'median_psf_fwhm', 'n_ccds', 'n_good_ccds', 'efftime']
exp_dr9 = join(exp_dr9, exp_all[exp_columns], keys='expnum', join_type='left')
print('exp_dr9', len(exp_dr9))

################################# Downselect the exposures #################################

np.random.seed(9371)
nexp_dict = {'g': 40, 'r': 35, 'z': 100}

expnum_list = []
print('Downselection')
for band in ['g', 'r', 'z']:
    idx = np.where(exp_dr9['filter']==band)[0]
    print(band, np.sum(exp_cosmos['filter']==band))
    idx = np.random.choice(idx, size=nexp_dict[band], replace=False)
    for index in idx:
        mask = exp_cosmos['filter']==band
        mask &= ~np.in1d(exp_cosmos['expnum'], expnum_list)
        index_best = np.argmin(np.log10(np.abs(exp_cosmos['efftime'][mask]/exp_dr9['efftime'][index])))
        expnum_list.append(exp_cosmos['expnum'][mask][index_best])

mask = np.in1d(exp_cosmos['expnum'], expnum_list)
exp_cosmos = exp_cosmos[mask]
print('exp_cosmos', len(exp_cosmos))

mask = (exp_cosmos['filter']!='r') | (exp_cosmos['efftime']<100)
exp_cosmos = exp_cosmos[mask]
print('exp_cosmos', len(exp_cosmos))

mask = np.in1d(ccd['expnum'], exp_cosmos['expnum'])
ccd = ccd[mask]
print('ccd', np.sum(mask))

# Flag the inner 1 deg CCDs
search_radius = 1.0*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(deep_ra, deep_dec, ccd['ra'], ccd['dec'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
ccd['inner'] = False
ccd['inner'][idx2] = True

################################## random subsets ##################################

seeds = {'g': 188, 'r': 538, 'z': 101}

idx_split = {}
efftime_list = {}
exp = {}
min_subset = 100000
for band in ['g', 'r', 'z']:

    np.random.seed(seeds[band])

    _, idx = np.unique(ccd['expnum'], return_index=True)
    exp[band] = ccd[idx].copy()

    mask = exp[band]['filter']==band
    exp[band] = exp[band][mask]

    # Shuffle
    idx = np.random.choice(len(exp[band]), size=len(exp[band]), replace=False)
    exp[band] = exp[band][idx]

    if band=='r':
        split_idx = np.cumsum(np.random.choice([1, 2, 3], p=[0.4, 0.5, 0.1], size=1000))
    else:
        split_idx = np.cumsum(np.random.choice([2, 3, 4, 5], p=[0.1, 0.3, 0.25, 0.35], size=1000))
    split_idx = split_idx[split_idx<len(idx)]
    if (split_idx[-1]!=len(idx)-1) and (split_idx[-1]<len(idx)-3):
        split_idx = np.append(split_idx, len(idx)-1)
    split = np.array_split(idx, split_idx)
    efftime_list[band] = np.zeros(len(split))
    for index in range(len(split)):
        mask = np.in1d(ccd['expnum'], exp[band]['expnum'][split[index]])
        mask &= ccd['inner']
        efftime_list[band][index] = np.sum(ccd['efftime'][mask])
    sort_subsets = np.argsort(efftime_list[band])
    efftime_list[band] = np.sort(efftime_list[band])
    idx_split[band] = [split[index] for index in sort_subsets]
    min_subset = np.minimum(min_subset, len(idx_split[band]))
    print('split', band, len(idx_split[band]))

print()
print('min_subset', min_subset)
print()

seeds = {'g': 188, 'r': 129, 'z': 116}

for band in ['g', 'r', 'z']:
    efftime_list[band] = list(efftime_list[band])
    np.random.seed(seeds[band])
    print(band, len(idx_split[band]))
    while len(idx_split[band])>min_subset:
        index_pop = np.random.choice(len(idx_split[band]))
        idx_split[band].pop(index_pop)
        efftime_list[band].pop(index_pop)
    print(len(idx_split[band]))

for index in range(min_subset):
    mask = np.full(len(ccd), False)
    for band in ['g', 'r', 'z']:
        mask |= np.in1d(ccd['expnum'], exp[band]['expnum'][idx_split[band][index]])
    ccd_subset = ccd[mask].copy()
    # ccd_subset.remove_columns('inner')
    ccd_subset.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/cosmos/survey-ccds-dr10-deep-fields-v1-defringed-subset-{}.fits'.format(index), overwrite=True)


# Sanity check: each CCD is only used once
fns = glob.glob('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/cosmos/*.fits')
print(len(fns))
cat = []
for fn in fns:
    tmp = Table(fitsio.read(fn))
    cat.append(tmp)
cat = vstack(cat)
cat['ccd_id_str'] = np.char.add(np.array(cat['expnum']).astype(str), cat['ccdname'])
if len(cat)!=len(np.unique(cat['ccd_id_str'])):
    raise ValueError

