from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

from astropy.coordinates import SkyCoord
from astropy import units as u


ext_coeffs = {
'g': 3.240, 'r': 2.276, 'i': 1.633, 'z': 1.263, 'y': 1.076,
'M411': 4.290, 'M438': 4.099, 'M464': 3.877, 'M490': 3.634, 'M517': 3.389}

# [decam-chatter 18711] Suggested MASKBITS for IBIS
# 2024-12-20
def get_maskbits_cut_bits(bands):

    cut_names = ['NPRIMARY', 'BRIGHT', 'MEDIUM', 'GALAXY']
    for band in bands:
        cut_names += ['SATUR_{}'.format(band), 'ALLMASK_{}'.format(band)]

    # pick one tractor file, to grab header
    fn = '/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis3/tractor/032/tractor-0329m042.fits'
    hdr = fitsio.read_header(fn, 0)
    keys = np.array([key for key in hdr if key[:5] == 'MBIT_'])
    names = np.array([hdr[key] for key in keys])

    cut_keys = keys[np.in1d(names, cut_names)]
    cut_bits = np.array([int(key.split('_')[-1]) for key in cut_keys])

    return cut_bits


field_ra, field_dec, field_rad_deg = 35.75, -4.75, 1.8

density_rad_deg = 1.6  # radius for density numbers
area = np.pi * density_rad_deg ** 2

##############################################################################################################################

ibis = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_xmm.fits'))
hsc = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_xmm-hsc_wide.fits'))
assert np.all(ibis['ibis_id']==hsc['ibis_id'])

print(len(ibis), len(hsc))
print(np.intersect1d(ibis.colnames, hsc.colnames))
hsc.remove_columns(['ibis_id', 'ra', 'dec'])
cat = hstack([ibis, hsc])
print(len(cat))

del ibis, hsc

# compute magnitudes
fnodet = 0.001
magnodet = 22.5 - 2.5 * np.log10(fnodet) # 30

# clip the magnitudes at 30.0
ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']
for band in ibis_filters:
    cat[band+'mag'] = 22.5 - 2.5*np.log10(np.clip(cat['flux_{}'.format(band)]*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))
    cat[band+'fibmag'] = 22.5 - 2.5*np.log10(np.clip(cat['fiberflux_{}'.format(band)]*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))
for band in ['g', 'r', 'i', 'z', 'y']:
    mask = ~np.isnan(cat['{}_cmodel_flux'.format(band)])
    cat[band+'mag'] = -99.
    cat[band+'mag'][mask] = 22.5 - 2.5*np.log10(np.clip(cat['{}_cmodel_flux'.format(band)]/3630.78*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))[mask]
    cat[band+'mag'][~mask] = 30.0

# leave-one-out average magnitude
for band in ibis_filters:
    cat[band+'_free_flux'] = 0.
    for band1 in ibis_filters:
        if band1!=band:
            cat[band+'_free_flux'] += cat['flux_'+band1] * 10**(0.4*ext_coeffs[band1]*cat['ebv']) / 4
    cat[band+'_free_mag'] = 22.5 - 2.5*np.log10(np.clip(cat[band+'_free_flux'], fnodet, None))

################### Add DESI LAE classifications ###################
desi = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/catalogs/tertiary49_ibis_lae_truth-20260127-add_target_info-add_lyaflux.fits'))
print(len(desi))

desi.rename_columns(['lae_sel', 'lae_sel_bright', 'M411_sel', 'M438_sel', 'M464_sel', 'M490_sel', 'M517_sel'], ['lae_sel_t49', 'lae_sel_bright_t49', 'M411_sel_t49', 'M438_sel_t49', 'M464_sel_t49', 'M490_sel_t49', 'M517_sel_t49'])

# https://github.com/rongpu/Python/blob/master/user_modules/match_coord.py
sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

idx1, idx2, d2d, d_ra, d_dec = match_coord(cat['ra'], cat['dec'], desi['RA'], desi['DEC'], search_radius=0.5, plot_q=True)
print(len(idx2), len(idx2)/len(desi))
cat['Z'] = -99.
cat['Z'][idx1] = desi['Z'][idx2]
cat['obs'] = False
cat['obs'][idx1] = True
cat['lae'] = False
cat['lae'][idx1] = desi['lae'][idx2]
cat['lae'] = False
cat['lae'][idx1] = desi['lae'][idx2]
for col in ['lae_sel_t49', 'lae_sel_bright_t49', 'M411_sel_t49', 'M438_sel_t49', 'M464_sel_t49', 'M490_sel_t49', 'M517_sel_t49']:
    cat[col] = False
    cat[col][idx1] = desi[col][idx2]
for col in ['DELTACHI2', 'lya_flux', 'lya_flux_err']:
    cat[col] = -99.
    cat[col][idx1] = desi[col][idx2]
print(len(desi), np.sum(cat['obs']))
print(np.sum(desi['lae']), np.sum(cat['lae']))

lyman_alpha_z = {'M411': [2.3, 2.51], 'M438': [2.51, 2.72], 'M464': [2.72, 2.93], 'M490': [2.93, 3.15], 'M517': [3.15, 3.36]}
lyman_alpha_z_min, lyman_alpha_z_max = 2.3, 3.36

cat['lae'] &= (cat['Z']>lyman_alpha_z_min) & (cat['Z']<lyman_alpha_z_max)
cat['lae_bright'] = np.full(len(cat), False)
for band in lyman_alpha_z.keys():
    zmin, zmax = lyman_alpha_z[band][0], lyman_alpha_z[band][1]
    mask = (cat['Z']>zmin) & (cat['Z']<=zmax)
    cat['lae_'+band] = cat['lae'] & mask
    cat['lae_'+band+'_bright'] = cat['lae'] & mask & (cat[band+'fibmag']<25.0)
    cat['lae_bright'] |= cat['lae_'+band+'_bright']

##############################################################################################################################

# cut on radius
field_c = SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg, frame='icrs')                                                                                                   
cs = SkyCoord(ra=cat['ra'] * u.deg, dec=cat['dec'] * u.deg, frame='icrs')
cat['CENTER_OFFSET'] = cs.separation(field_c).deg
mask = cat['CENTER_OFFSET'] < field_rad_deg
cat = cat[mask]
print('radius cut', len(cat), np.sum(mask)/len(mask))

mask_rad = cat['CENTER_OFFSET'] < density_rad_deg
print('{:.1f}% within {:.2f} deg radius'.format(100*np.sum(mask_rad)/len(mask_rad), density_rad_deg))

##############################################################################################################################

def select_lae(cat, band, cfg, area):
    print(f'\n{band}')
    mask = np.full(len(cat), True)
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tinitial')

    # cut on maskbits
    cut_bits = get_maskbits_cut_bits(ibis_filters)
    for bit in cut_bits:
        mask &= (cat['maskbits'] & 2**bit) == 0
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tafter maskbits cut')

    # cut on IBIS fibmag
    fibmag = cat[f'{band}fibmag']
    mask &= (fibmag > FIBMAG_MIN) & (fibmag < FIBMAG_MAX)
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tafter {FIBMAG_MIN} < fibmag < {FIBMAG_MAX} cut')

    # cut on IBIS color
    color = cat[band+'_free_mag'] - cat[band+'mag']
    mask &= color > np.exp((cat[f'{band}fibmag'] - 25.) * cfg['exp_a']) * cfg['exp_b'] + cfg['mb_colmin']
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tafter IBIS color cut')

    # cut on g-r and r-i
    mask1 = (cat['gmag'] - cat['rmag'] < cfg['gr_colmax']) | ((cat['gmag'] < MAGNODET) & (cat['rmag'] == MAGNODET)) | ((cat['gmag'] == MAGNODET) & (cat['rmag'] == MAGNODET))
    print(f'{(mask & mask1 & mask_rad).sum() / area:.1f}/deg2\tafter g-r < {cfg["gr_colmax"]} cut')
    mask1 &= (cat['rmag'] - cat['imag'] < cfg['ri_colmax']) | ((cat['rmag'] < MAGNODET) & (cat['imag'] == MAGNODET)) | ((cat['rmag'] == MAGNODET) & (cat['imag'] == MAGNODET))
    print(f'{(mask & mask1 & mask_rad).sum() / area:.1f}/deg2\tafter r-i < {cfg["ri_colmax"]} cut')
    mask1 |= (cat['rmag'] > MAGNODET_HSC)
    print(f'{(mask & mask1 & mask_rad).sum() / area:.1f}/deg2\tafter adding back rmag>{MAGNODET_HSC} sources')
    mask &= mask1

    # remove shreds of largish galaxies and brightish stars
    mask &= cat['fracflux_'+band]<1
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tafter removing shreds')

    # better safe than sorry... remove a handful of dubious bright r-mag targets
    mask1 = cat['rmag'] < RMAG_BRIGHT_CUT
    print(f'remove {np.sum(mask & mask1)} targets with r_mag < {RMAG_BRIGHT_CUT}')
    mask &= (~mask1)

    print('LAE completeness: {:.1f}% ({}/{})'.format(100*np.sum(cat['lae_'+band+'_bright'][mask])/np.sum(cat['lae_'+band+'_bright']), np.sum(cat['lae_'+band+'_bright'][mask]), np.sum(cat['lae_'+band+'_bright'])))
    print('LAE purity:       {:.1f}% ({}/{})'.format(100*np.sum(cat['lae'][mask])/np.sum(cat['obs'][mask]), np.sum(cat['lae'][mask]), np.sum(cat['obs'][mask])))

    return mask


FIBMAG_MIN, FIBMAG_MAX = 22.0, 25.0
RMAG_BRIGHT_CUT = 18
MAGNODET = 30.0
MAGNODET_HSC = 25.0

LAE_CONFIGS = {
    'M411': {
        'mb_colmin': 0.24,
        'exp_a': 4.,
        'exp_b': 0.3,
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
    },
    'M438': {
        'mb_colmin': 0.23,
        'exp_a': 4.,
        'exp_b': 0.6,
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
    },
    'M464': {
        'mb_colmin': 0.30,
        'exp_a': 4.,
        'exp_b': 0.3,
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
    },
    'M490': {
        'mb_colmin': 0.55,
        'exp_a': 4.,
        'exp_b': 0.3,
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
    },
    'M517': {
        'mb_colmin': 0.81,
        'exp_a': 4.,
        'exp_b': 0.3,
        'ri_colmax': 0.4, 'gr_colmax': 0.8,
    },
}

for band, cfg in LAE_CONFIGS.items():
    cat[f'{band}_sel'] = select_lae(cat, band, cfg, area)

print('\nAll combined:')
cat['lae_sel'] = False
for band in ibis_filters:
    cat['lae_sel'] |= cat[f'{band}_sel']
    print(f'{np.sum(cat["lae_sel"] & mask_rad) / area:.1f}/deg2\tafter adding {band}')
mask = cat['lae_sel'].copy()
print('LAE completeness: {:.1f}% ({}/{})'.format(100*np.sum(cat['lae_bright'][mask])/np.sum(cat['lae_bright']), np.sum(cat['lae_bright'][mask]), np.sum(cat['lae_bright'])))
print('LAE purity:       {:.1f}% ({}/{})'.format(100*np.sum(cat['lae'][mask])/np.sum(cat['obs'][mask]), np.sum(cat['lae'][mask]), np.sum(cat['obs'][mask])))

# cat.write('/pscratch/sd/r/rongpu/tmp/ibis_hsc-wide_xmm-synthmag.fits', overwrite=False)
cat.write('/pscratch/sd/r/rongpu/tmp/tertiary54/tertiary54_lae_targets-xmm-all.fits', overwrite=False)

mask = cat['lae_sel'].copy()
cat[mask].write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/tertiary54_lae_targets-xmm.fits', overwrite=False)
