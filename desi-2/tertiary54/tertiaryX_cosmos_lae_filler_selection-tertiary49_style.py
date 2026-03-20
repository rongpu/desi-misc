# Running the same selection (except for the IBIS color cuts which are adjusted to conserve the target density) in COSMOS; for testing purposes only

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

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
    fn = '/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis3/tractor/032/tractor-0329m042.fits'
    hdr = fitsio.read_header(fn, 0)
    keys = np.array([key for key in hdr if key[:5] == 'MBIT_'])
    names = np.array([hdr[key] for key in keys])

    cut_keys = keys[np.in1d(names, cut_names)]
    cut_bits = np.array([int(key.split('_')[-1]) for key in cut_keys])

    return cut_bits


# offset field by +0.2 in dec to avoid shallow m464 region
field_ra, field_dec, field_rad_deg = 150.1, 2.182 + 0.2, 1.8
area = np.pi * field_rad_deg ** 2

################################################################################################################################################
# Load IBIS tractor catalog and cross-matche HSC-deep and HSC-wide

ibis = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_cosmos.fits'))
hsc_deep = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_cosmos-hsc_deep.fits'))
hsc_wide = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/tractor_cosmos-hsc_wide.fits'))
assert np.all(ibis['ibis_id']==hsc_deep['ibis_id']) and np.all(ibis['ibis_id']==hsc_wide['ibis_id'])

cols = list(np.intersect1d(hsc_deep.colnames, hsc_wide.colnames))
hsc_deep = hsc_deep[cols]
hsc_wide = hsc_wide[cols]
hsc_deep['hsc_deep'] = True
hsc_wide['hsc_wide'] = False
hsc_deep['id'] = np.arange(len(hsc_deep))
hsc_wide['id'] = np.arange(len(hsc_wide))

# Use HSC-deep if available
mask_deep = hsc_deep['ra']!=-99
hsc = vstack([hsc_deep[mask_deep], hsc_wide[~mask_deep]])
hsc.sort('id')
hsc.remove_column('id')
assert len(hsc)==len(hsc_deep) and np.all(hsc['ibis_id']==hsc_deep['ibis_id'])

print(len(ibis), len(hsc))
print(np.intersect1d(ibis.colnames, hsc.colnames))
hsc.remove_columns(['ibis_id', 'ra', 'dec'])
cat = hstack([ibis, hsc])
print(len(cat))

################################################################################################################################################

# cut on radius
field_c = SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg, frame='icrs')                                                                                                   
cs = SkyCoord(ra=cat['ra'] * u.deg, dec=cat['dec'] * u.deg, frame='icrs')
cat['FIELD_RADIUS_OFFSET'] = cs.separation(field_c).deg
mask = cat['FIELD_RADIUS_OFFSET'] < field_rad_deg
print('{:.1f}k/deg2\tafter radius cut'.format(mask.sum() / 1000. / area))
cat = cat[mask]
print(len(cat), np.sum(mask)/len(mask))

################################################################################################################################################

# compute magnitudes
fnodet = 0.001
magnodet = 22.5 - 2.5 * np.log10(fnodet) # 30

# clip the magnitudes at 30.0
for band in ['M411', 'M438', 'M464', 'M490', 'M517']:
    cat[band+'fibmag'] = 22.5 - 2.5*np.log10(np.clip(cat['fiberflux_{}'.format(band)]*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))
for band in ['g', 'r', 'i', 'z', 'y']:
    mask = ~np.isnan(cat['{}_cmodel_flux'.format(band)])
    cat[band+'mag'] = -99.
    cat[band+'mag'][mask] = 22.5 - 2.5*np.log10(np.clip(cat['{}_cmodel_flux'.format(band)]/3630.78*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))[mask]
    cat[band+'mag'][~mask] = 30.0

##################### LAE selections #####################

LAE_CONFIGS = {
    'M411': {
        'neighbors': ('M411', 'M438'),
        'blue_band': 'M438', 'red_band': 'M411',
        'mb_colmin': 0.1,
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
    },
    'M438': {
        'neighbors': ('M411', 'M438'),
        'blue_band': 'M411', 'red_band': 'M438',
        'mb_colmin': 0.31,
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
    },
    'M464': {
        'neighbors': ('M438', 'M464'),
        'blue_band': 'M438', 'red_band': 'M464',
        'mb_colmin': 0.3,
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
    },
    'M490': {
        'neighbors': ('M464', 'M490'),
        'blue_band': 'M464', 'red_band': 'M490',
        'mb_colmin': 0.32,
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
    },
    'M517': {
        'neighbors': ('M490', 'M517'),
        'blue_band': 'M490', 'red_band': 'M517',
        'mb_colmin': 0.29,
        'ri_colmax': 0.4,  'gr_colmax': 0.8,
    },
}

FIBMAG_MIN, FIBMAG_MAX = 22.0, 25.25
FIBMAG_MAX_BRIGHT = 25.0
RMAG_BRIGHT_CUT = 18
MAGNODET_HSC = 25.0


def select_lae(cat, band, cfg, area, magnodet):
    print(f'\n{band}')
    mask = np.full(len(cat), True)

    # cut on maskbits
    cut_bits = get_maskbits_cut_bits(list(cfg['neighbors']))
    for bit in cut_bits:
        mask &= (cat['maskbits'] & 2**bit) == 0
    print(f'{mask.sum() / 1000. / area:.1f}k/deg2\tafter maskbits cut')

    # cut on IBIS fibmag
    fibmag = cat[f'{band}fibmag']
    mask &= (fibmag > FIBMAG_MIN) & (fibmag < FIBMAG_MAX)
    print(f'{mask.sum() / 1000. / area:.1f}k/deg2\tafter {FIBMAG_MIN} < fibmag < {FIBMAG_MAX} cut')

    # cut on IBIS color
    color = cat[f'{cfg["blue_band"]}fibmag'] - cat[f'{cfg["red_band"]}fibmag']
    mask &= color > cfg['mb_colmin']
    print(f'{mask.sum() / 1000. / area:.1f}k/deg2\tafter IBIS color > {cfg["mb_colmin"]} cut')

    # cut on g-r and r-i
    mask1 = (cat['gmag'] - cat['rmag'] < cfg['gr_colmax']) | ((cat['gmag'] < magnodet) & (cat['rmag'] == magnodet))
    print(f'{(mask & mask1).sum() / 1000. / area:.1f}k/deg2\tafter g-r < {cfg["gr_colmax"]} cut')
    mask1 &= (cat['rmag'] - cat['imag'] < cfg['ri_colmax']) | ((cat['rmag'] < magnodet) & (cat['imag'] == magnodet))
    print(f'{(mask & mask1).sum() / 1000. / area:.1f}k/deg2\tafter r-i < {cfg["ri_colmax"]} cut')
    mask1 |= (cat['gmag'] > MAGNODET_HSC) | (cat['rmag'] > MAGNODET_HSC) | (cat['imag'] > MAGNODET_HSC)  # also keep everything fainter than mag=25
    mask &= mask1
    print(f'{mask.sum() / 1000. / area:.1f}k/deg2\tafter adding back faint sources')

    # better safe than sorry... remove a handful of dubious bright r-mag targets
    mask1 = cat['rmag'] < RMAG_BRIGHT_CUT
    print(f'remove {np.sum(mask & mask1)} targets with r_mag < {RMAG_BRIGHT_CUT}')
    mask &= (~mask1)

    # brighter IBIS limit
    mask_bright = mask.copy()
    mask_bright &= cat[f'{band}fibmag'] < FIBMAG_MAX_BRIGHT
    print('{:.1f}k/deg2\tafter fibmag_brightcut'.format(mask_bright.sum() / 1000. / area))

    return mask, mask_bright


for band, cfg in LAE_CONFIGS.items():
    cat[f'{band}_sel'], cat[f'{band}_sel_bright'] = select_lae(cat, band, cfg, area, magnodet)

print('\nAll combined:')
cat['lae_sel'] = False
for band in LAE_CONFIGS:
    cat['lae_sel'] |= cat[f'{band}_sel']
    print(f'{np.sum(cat["lae_sel"]) / 1000. / area:.1f}k/deg2\tafter adding {band}')

print('\nAll bright combined:')
cat['lae_sel_bright'] = False
for band in LAE_CONFIGS:
    cat['lae_sel_bright'] |= cat[f'{band}_sel_bright']
    print(f'{np.sum(cat["lae_sel_bright"]) / 1000. / area:.1f}k/deg2\tafter adding {band}')

assert np.sum(cat['lae_sel_bright'])==np.sum(cat['lae_sel_bright'] & cat['lae_sel'])
print()

mask = cat['lae_sel'].copy()
cat = cat[mask]
print(len(cat))

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/tertiary54_cosmos_lae_filler_targets.fits', overwrite=False)
