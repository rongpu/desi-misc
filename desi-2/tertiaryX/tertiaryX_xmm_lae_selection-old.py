# High priority LAEs selected using IBIS-DR1 (deep imaging) and HSC forced photometry
# Note the change of IBIS extinction coefficients

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
from scipy.stats import norm


ext_coeffs = {
'g': 3.240, 'r': 2.276, 'i': 1.633, 'z': 1.263, 'y': 1.076,
'M411': 3.774, 'M438': 3.576, 'M464': 3.366, 'M490': 3.141, 'M517': 2.922}

# [decam-chatter 18711] Suggested MASKBITS for IBIS
# 2024-12-20
def get_maskbits_cut_bits(bands):

    cut_names = ['NPRIMARY', 'BRIGHT', 'MEDIUM', 'GALAXY']
    for band in bands:
        cut_names += ['SATUR_{}'.format(band), 'ALLMASK_{}'.format(band)]

    # pick one tractor file, to grab header
    fn = '/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/tractor/032/tractor-0329m042.fits'
    hdr = fitsio.read_header(fn, 0)
    keys = np.array([key for key in hdr if key[:5] == 'MBIT_'])
    names = np.array([hdr[key] for key in keys])

    cut_keys = keys[np.in1d(names, cut_names)]
    cut_bits = np.array([int(key.split('_')[-1]) for key in cut_keys])

    return cut_bits


def flux_threshold(f_b, depth, c0, prob_tail, rescale_factor=1, include_color_lim=False):
    '''
    Compute the continuum flux threshold f_free_lim for selecting narrow-band excess objects (e.g. LAEs).

    This threshold is used to select objects with f_free < f_free_lim, where f_free is the measured
    continuum flux: the mean of the 4 other IBIS bands.

    Contaminant model: every contaminant has the same true color
    c0 = free_mag - mag, i.e. true band flux = a * true continuum flux with
    a = 10**(0.4*c0). For a contaminant, the excess X = f_b - a*f_free then has
    mean zero and Gaussian error sigma_x, whatever its brightness. Requiring
    X > z*sigma_x (z = norm.isf(prob_tail)) therefore lets each contaminant pass
    with probability prob_tail at every magnitude.

    Assumes Gaussian, background-limited flux errors that are independent
    between bands, and the same depth in all bands.

    Parameters
    ----------
    f_b : float or array
        Measured selection-band flux (total, extinction-corrected; nanomaggies).
    depth : float
        5-sigma point-source depth (mag), the same for all bands.
    c0 : float
        True color (free_mag - mag) of the contaminants.
    prob_tail : float
        Probability that a contaminant passes the cut.
    rescale_factor : float
        Multiplier on the flux errors of all bands.
    include_color_lim : bool
        If True, also return the equivalent color threshold.

    Returns
    -------
    f_free_lim : array
        Upper limit on f_free (nanomaggies). Compare against the unclipped
        continuum flux, not one derived from a clipped free_mag.
    color_lim : array, only if include_color_lim is True
        Equivalent lower limit on color = free_mag - mag; select
        color > color_lim. It is inf where f_free_lim <= 0 < f_b (no object
        with positive f_free passes) and NaN where f_b <= 0 (color undefined).
    '''
    f_b = np.asarray(f_b, dtype=float)
    z = norm.isf(prob_tail)  # gaussian sigma at tail probability prob_tail
    a = 10**(0.4 * c0)
    sigma_b = rescale_factor * 10**(-0.4 * (depth - 22.5)) / 5  # error on selection-band flux
    sigma_x = sigma_b * np.sqrt(1 + a**2 / 4)  # error on f_b - a*f_free
    f_free_lim = (f_b - z * sigma_x) / a

    if not include_color_lim:
        return f_free_lim
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            color_lim = np.where(f_free_lim > 0, 2.5 * np.log10(f_b / f_free_lim), np.inf)
        color_lim = np.where(f_b > 0, color_lim, np.nan)
        return f_free_lim, color_lim


field_ra, field_dec, field_rad_deg = 35.75, -4.75, 1.8

density_rad_deg = 1.6  # radius for density numbers
area = np.pi * density_rad_deg ** 2

##############################################################################################################################

ibis = Table(fitsio.read('/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-xmm-combined.fits'))
hsc = Table(fitsio.read('/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-xmm-hsc-wide-forced.fits'))
assert np.all(ibis['ibis_id_dr1']==hsc['ibis_id_dr1'])

print(len(ibis), len(hsc))
print(np.intersect1d(ibis.colnames, hsc.colnames))
hsc.remove_columns(['ibis_id_dr1'])
cat = hstack([ibis, hsc])
print(len(cat))

del ibis, hsc

# cut on radius
field_c = SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg, frame='icrs')
cs = SkyCoord(ra=cat['ra'] * u.deg, dec=cat['dec'] * u.deg, frame='icrs')
cat['CENTER_OFFSET'] = cs.separation(field_c).deg
mask = cat['CENTER_OFFSET'] < field_rad_deg
cat = cat[mask]
print('radius cut', len(cat), np.sum(mask)/len(mask))

##############################################################################################################################
# compute magnitudes

# clip the magnitudes at 30.0
fnodet = 0.001
magnodet = 22.5 - 2.5 * np.log10(fnodet)  # 30

ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']
for band in ibis_filters:
    cat[band+'mag'] = 22.5 - 2.5*np.log10(np.clip(cat['flux_{}'.format(band)]*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))
    cat[band+'fibmag'] = 22.5 - 2.5*np.log10(np.clip(cat['fiberflux_{}'.format(band)]*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))

for band in ['g', 'r', 'i', 'z', 'y']:
    cat[band+'mag'] = 22.5 - 2.5*np.log10(np.clip(cat['forced_flux_{}'.format(band)]*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))

# leave-one-out average magnitude
for band in ibis_filters:
    cat[band+'_free_flux'] = 0.
    for band1 in ibis_filters:
        if band1!=band:
            cat[band+'_free_flux'] += cat['flux_'+band1] * 10**(0.4*ext_coeffs[band1]*cat['ebv']) / 4
    cat[band+'_free_mag'] = 22.5 - 2.5*np.log10(np.clip(cat[band+'_free_flux'], fnodet, None))

# compute depths
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    for band in ibis_filters:
        cat[band+'_depth_ebv'] = -2.5*(np.log10((5/np.sqrt(cat['psfdepth_'+band])))-9) - ext_coeffs[band]*cat['ebv']

##############################################################################################################################

mask_rad = cat['CENTER_OFFSET'] < density_rad_deg
print('{:.1f}% within {:.2f} deg radius'.format(100*np.sum(mask_rad)/len(mask_rad), density_rad_deg))

median_depths = {}
print('Median depths:')
for band in ibis_filters:
    median_depths[band] = np.nanmedian(cat[band+'_depth_ebv'][mask_rad])
    print(band, '{:.2f}'.format(median_depths[band]))

cut_bits = get_maskbits_cut_bits(ibis_filters)

##############################################################################################################################


def select_lae(cat, band, cfg, area):
    print(f'\n{band}')
    mask = np.full(len(cat), True)
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tinitial')

    # cut on maskbits
    for bit in cut_bits:
        mask &= (cat['maskbits'] & 2**bit) == 0
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tafter maskbits cut')

    # cut on IBIS fibmag
    fibmag = cat[f'{band}fibmag']
    mask &= (fibmag > FIBMAG_MIN) & (fibmag < FIBMAG_MAX)
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tafter {FIBMAG_MIN} < fibmag < {FIBMAG_MAX} cut')

    # cut on g-r and r-i
    mask1 = (cat['gmag'] - cat['rmag'] < cfg['gr_colmax']) | ((cat['gmag'] < MAGNODET) & (cat['rmag'] == MAGNODET)) | ((cat['gmag'] == MAGNODET) & (cat['rmag'] == MAGNODET))
    print(f'{(mask & mask1 & mask_rad).sum() / area:.1f}/deg2\tafter g-r < {cfg["gr_colmax"]} cut')
    mask1 &= (cat['rmag'] - cat['imag'] < cfg["ri_colmax"]) | ((cat['rmag'] < MAGNODET) & (cat['imag'] == MAGNODET)) | ((cat['rmag'] == MAGNODET) & (cat['imag'] == MAGNODET))
    print(f'{(mask & mask1 & mask_rad).sum() / area:.1f}/deg2\tafter r-i < {cfg["ri_colmax"]} cut')
    mask1 |= (cat['rmag'] > MAGNODET_HSC)
    print(f'{(mask & mask1 & mask_rad).sum() / area:.1f}/deg2\tafter adding back rmag>{MAGNODET_HSC} sources')
    mask &= mask1

    # remove shreds of largish galaxies and brightish stars
    mask &= cat['fracflux_'+band]<1
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tafter removing shreds')

    # cut on excess IBIS flux in the selection band
    color_offset = cfg['color_offset']
    median_depth = median_depths[band]
    mask1 = mask.copy() & (cat[f'{band}mag']>24.4) & (cat[f'{band}mag']<24.5) & mask_rad
    median_color = np.nanmedian((cat[band+'_free_mag'] - cat[band+'mag'])[mask1])
    print('median IBIS color = {:.2f}'.format(median_color))
    flux = cat['flux_{}'.format(band)]*10**(0.4*ext_coeffs[band]*cat['ebv'])
    mask &= cat[band+'_free_flux'] < 10**(-color_offset/2.5) * flux_threshold(flux, depth=median_depth, c0=median_color, prob_tail=0.0001, include_color_lim=False)
    mask &= cat[band+'_free_mag'] - cat[band+'mag'] > cfg['mb_colmin']
    print(f'{(mask & mask_rad).sum() / area:.1f}/deg2\tafter IBIS color cut')

    # better safe than sorry... remove a handful of dubious bright r-mag targets
    mask1 = cat['rmag'] < RMAG_BRIGHT_CUT
    print(f'remove {np.sum(mask & mask1)} targets with r_mag < {RMAG_BRIGHT_CUT}')
    mask &= (~mask1)

    return mask


FIBMAG_MIN, FIBMAG_MAX = 22.0, 25.0
RMAG_BRIGHT_CUT = 18
MAGNODET = 30.0
MAGNODET_HSC = 25.0

LAE_CONFIGS = {
    'M411': {
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
        'color_offset': 0.07,
        'mb_colmin': 0.15,
    },
    'M438': {
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
        'color_offset': 0.06,
        'mb_colmin': 0.15,
    },
    'M464': {
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
        'color_offset': 0.08,
        'mb_colmin': 0.25,
    },
    'M490': {
        'ri_colmax': 0.35, 'gr_colmax': 0.7,
        'color_offset': 0.21,
        'mb_colmin': 0.45,
    },
    'M517': {
        'ri_colmax': 0.4, 'gr_colmax': 0.8,
        'color_offset': 0.45,
        'mb_colmin': 0.7,
    },
}

for band, cfg in LAE_CONFIGS.items():
    cat[f'{band}_sel'] = select_lae(cat, band, cfg, area)

print('\nAll combined:')
cat['lae_sel'] = False
for band in ibis_filters:
    cat['lae_sel'] |= cat[f'{band}_sel']
    band_density = np.sum(cat[f'{band}_sel'] & mask_rad) / area
    print('{:8.1f}/deg2 after adding {} ({:.1f}/deg2)'.format(np.sum(cat['lae_sel'] & mask_rad) / area, band, band_density))


cat.write('/pscratch/sd/r/rongpu/tmp/tertiaryX/tertiaryX_lae_targets-all.fits', overwrite=True)

mask = cat['lae_sel'].copy()
cat[mask].write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_lae_targets.fits', overwrite=False)
