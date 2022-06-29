# CCDs with ra, dec within 2.2 deg of the field centers

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord, match_self


# bad exposures found in visual inspection
bad_expids = [938452, 964237, 964238, 964239, 980290, 980291, 986155, 1002854, 1006183, 1006184]

###################### Exposure-level cuts ######################

exp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures-medians.fits'))
print(len(exp))

columns = ['image_filename', 'camera', 'expnum', 'plver', 'procdate', 'plprocid', 'object', 'propid', 'filter', 'exptime', 'mjd_obs', 'airmass', 'ra_bore', 'dec_bore', 'zpt']
tmp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-unique-exposures.fits', columns=columns))
print(len(tmp))

exp = join(exp, tmp, keys='expnum')

mask = exp['filter']=='z'
exp = exp[mask]
print(len(exp))

mask = exp['n_good_ccds']==61
print(np.sum(mask)/len(mask))
exp = exp[mask]
print(len(exp))

mask = np.in1d(exp['expnum'], bad_expids)
exp = exp[~mask]
print(len(exp))

# Earlier PLVERs have uggly sawtooth patterns in the outer CCDs
mask = (exp['plver']>='V4.9')
exp = exp[mask]
print(len(exp))

# # Remove exposures with artifacts
# mask = np.char.startswith(exp['image_filename'], 'decam/CP/V4.8.2a/CP20131128')  # many of these exposures look bad
# exp = exp[~mask]
# print(len(exp))

# # select exposures after S30 came back to life
# mask = exp['mjd_obs']>57500
# exp = exp[mask]
# print(len(exp))

mask = exp['exptime']>=60
exp = exp[mask]
print(len(exp))

psfex = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/survey-ccds-dr10-v4-psfex-fwhm.fits'))
_, idx = np.unique(psfex['expnum'], return_index=True)
psfex_exp = psfex[idx].copy()
exp = join(exp, psfex_exp, join_type='inner', keys='expnum')
print(len(exp))

# Quality cuts for deep exposures
seeing_limits = {'g':1.65, 'r':1.55, 'i':1.4, 'z':1.4, 'Y':2.0}
ccdskycounts_limits = {'z':35}
zpt_limits = {'g':24.9, 'r':25.15, 'i':25.15, 'z':24.85, 'Y':23.7}

exp['quality'] = False
for band in ['z']:
    print(band)
    mask = exp['filter']==band
    print(np.sum(mask))
    exp['quality'][mask] = exp['median_psf_fwhm'][mask]*0.262<seeing_limits[band]
    print(np.sum(exp['quality'][mask]))
    exp['quality'][mask] &= exp['median_ccdskycounts'][mask]<ccdskycounts_limits[band]
    print(np.sum(exp['quality'][mask]))
    exp['quality'][mask] &= exp['zpt'][mask]>zpt_limits[band]
    print(np.sum(exp['quality'][mask]))

mask = exp['quality'].copy()
exp = exp[mask]
print('exp', len(exp))

# Impose dithering: no neighboring exposures within 600 arcsec
n_duplicates = np.inf
while n_duplicates>0:
    n_duplicates, idx1, idx2 = match_self(exp['ra_bore'], exp['dec_bore'], search_radius=600, return_indices=True, plot_q=False)
    mask_keep = np.full(len(exp), True)
    mask = exp['exptime'][idx1]<exp['exptime'][idx2]
    mask_keep[idx1[mask]] = False
    mask_keep[idx2[~mask]] = False
    exp = exp[mask_keep]
print('exp', len(exp))

# Remove exposures within 2 degrees of a mag<2 star
stars = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_reference_dr9.fits', columns=['mask_mag']))
idx = np.where(stars['mask_mag']<2)[0]
stars = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_reference_dr9.fits', rows=idx))
print('Bright stars', len(stars))
search_radius = 2.*3600.  # arcsec
idx1, idx2, d2d, d_ra, d_dec = match_coord(stars['RA'], stars['DEC'], exp['ra_bore'], exp['dec_bore'], search_radius=search_radius, plot_q=False, keep_all_pairs=True)
mask = np.in1d(np.arange(len(exp)), idx2)
exp = exp[~mask]
print('exp', len(exp))

# # Require DR9 fringe for normalization
# dr9_fringe = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-z-fringe-stacking-dr9-status.fits'))
# exp = join(exp, dr9_fringe[['expnum', 'old_fringe', 'new_fringe']], join_type='inner')
# print('exp', len(exp))
# mask = exp['new_fringe'].copy()
# exp = exp[mask]
# print('exp', len(exp))

# Require DR10 blobmask
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_ccd_blob_mask'
mask = np.full(len(exp), True)
for index in range(len(exp)):
    # if index%100==0:
    #     print(index, len(exp))
    str_loc = str.find(exp['image_filename'][index].strip(), '.fits')
    img_filename_base = exp['image_filename'][index].strip()[:str_loc]
    blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
    if (not os.path.isfile(blob_path)) or os.stat(blob_path).st_size==0:
        # print(blob_path)
        mask[index] = False
exp = exp[mask]
print('exp with blob', len(exp))

###################### CCD-level cuts ######################

ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', columns=['expnum']))
mask = np.in1d(ccd['expnum'], exp['expnum'])
idx = np.where(mask)[0]

ccd = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-dr10-v4.fits', rows=idx))
print('ccd', len(ccd))

# mask = np.in1d(ccd['expnum'], exp['expnum'])
# ccd = ccd[mask]
# print('ccd', np.sum(mask))

# basic quality cuts
basic_quality = (ccd['ccd_cuts']==0) | (ccd['ccd_cuts']==2**14)  # ignore DEPTH_CUT
ccd = ccd[basic_quality]
print('ccd quality', len(ccd))

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

tmp = Table()
tmp['expnum'], tmp['count'] = np.unique(ccd['expnum'], return_counts=True)
mask = tmp['count']==61
print(np.sum(mask)/len(mask))
tmp = tmp[mask]
mask = np.in1d(ccd['expnum'], tmp['expnum'])
ccd = ccd[mask]
print(len(ccd))

###################### Get final CCD list ######################

exp = exp[['expnum', 'median_fwhm', 'median_ccdskycounts', 'median_psf_fwhm', 'n_ccds', 'n_good_ccds']]
ccd = join(ccd, exp, keys='expnum', join_type='left')
print('ccd', len(ccd))

ccd['efftime'] = 10**(0.4*ccd['ccdzpt']-9) * ccd['exptime'] / (ccd['median_ccdskycounts'] * ccd['psf_fwhm']**2)
# ccd['exposure_efftime'] = 10**(0.4*ccd['zpt']-9) * ccd['exptime'] / (ccd['median_ccdskycounts'] * ccd['median_psf_fwhm']**2)

ccd.rename_column('ccd_cuts', 'ccd_cuts_dr10')
ccd['ccd_cuts'] = 0

ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/survey-ccds-z-fringe.fits', overwrite=True)
