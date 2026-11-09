# High priority LAEs selected using IBIS-DR1 (deep imaging) and HSC forced photometry
# Note the change of IBIS extinction coefficients

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

# https://github.com/rongpu/Python/blob/master/user_modules/match_coord.py
sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

##############################################################################################################################
# Configuration

ibis_dir = '/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1'
ibis_fn = f'{ibis_dir}/catalogs/tractor-xmm-combined.fits'
hsc_fn = f'{ibis_dir}/catalogs/tractor-xmm-hsc-wide-forced.fits'
maskbits_header_fn = f'{ibis_dir}/tractor/032/tractor-0329m042.fits'  # any tractor file, to grab the MBIT_ header
output_all_fn = '/pscratch/sd/r/rongpu/tmp/tertiaryX/tertiaryX_lae_targets-all.fits'
output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-ibis_deep.fits'

ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']
hsc_filters = ['g', 'r', 'i', 'z', 'y']

ext_coeffs = {
    'g': 3.240, 'r': 2.276, 'i': 1.633, 'z': 1.263, 'y': 1.076,
    'M411': 3.774, 'M438': 3.576, 'M464': 3.366, 'M490': 3.141, 'M517': 2.922}

# From compute_contaminant_median_colors.py
contam_median_colors = {'M411': -0.07, 'M438': -0.02, 'M464': 0.03, 'M490': 0.1, 'M517': 0.12}

field_ra, field_dec, field_rad_deg = 35.75, -4.75, 1.8
density_rad_deg = 1.6  # radius for density numbers
area = np.pi * density_rad_deg ** 2

# clip the magnitudes at 30.0
fnodet = 0.001
magnodet = 22.5 - 2.5 * np.log10(fnodet)  # 30

FIBMAG_MIN, FIBMAG_MAX = 22.0, 25.0
RMAG_ADD_BACK_LIMIT = 25.0
RMAG_BRIGHT_CUT = 18
PROB_TAIL = 0.0001

LAE_CONFIGS = {
    'M411': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': -0.075, 'color_offset_sigma': 0.087},
    'M438': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': -0.085, 'color_offset_sigma': 0.08},
    'M464': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': -0.125, 'color_offset_sigma': 0.098},
    'M490': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': -0.03, 'color_offset_sigma': 0.121},
    'M517': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': 0., 'color_offset_sigma': 0.178},
}

##############################################################################################################################
# Functions


# [decam-chatter 18711] Suggested MASKBITS for IBIS 2024-12-20
def get_maskbits_cut_bits(bands, fn):
    cut_names = ['NPRIMARY', 'BRIGHT', 'MEDIUM', 'GALAXY']
    cut_names += [f'{prefix}_{band}' for band in bands for prefix in ['SATUR', 'ALLMASK']]
    print('MASKBITS names:', cut_names)

    hdr = fitsio.read_header(fn, 0)  # MBIT_<bit> = '<name>'
    name_to_bit = {hdr[key]: int(key[len('MBIT_'):]) for key in hdr if key.startswith('MBIT_')}

    missing = [name for name in cut_names if name not in name_to_bit]
    if missing:
        raise ValueError(f'maskbit names not found in {fn}: {missing}')

    return np.array(sorted(name_to_bit[name] for name in cut_names))


def flux_threshold(f_b, depth, c0, prob_tail, sigma_c=0., rescale_factor=1, include_color_lim=False):
    '''
    Compute the continuum flux threshold f_ref_lim for selecting narrow-band excess objects (e.g. LAEs).

    This threshold is used to select objects with f_ref < f_ref_lim, where f_ref is the measured
    continuum flux: the mean of the 4 other IBIS bands.

    Contaminant model: every contaminant has the same true color
    c0 = ref_mag - mag, i.e. true band flux = a * true continuum flux with
    a = 10**(0.4*c0). For a contaminant, the excess X = f_b - a*f_ref then has
    mean zero and Gaussian error sigma_x, whatever its brightness. Requiring
    X > z*sigma_x (z = norm.isf(prob_tail)) therefore lets each contaminant pass
    with probability prob_tail at every magnitude.

    Assumes Gaussian, background-limited flux errors that are independent
    between bands, and the same depth in all bands.

    Intrinsic color scatter: if sigma_c > 0, the contaminant true color is
    c ~ N(c0, sigma_c), which adds a flux-proportional term to the error on X:
    Var(X) = sigma_x**2 + (k * a * F)**2, where F is the true continuum flux
    (estimated as f_ref at the threshold). k is chosen so that the bright-end
    limit is exactly color > c0 + z*sigma_c. sigma_c = 0 recovers the
    single-color model.

    Parameters
    ----------
    f_b : float or array
        Measured selection-band flux (total, extinction-corrected; nanomaggies).
    depth : float
        5-sigma point-source depth (mag), the same for all bands.
    c0 : float
        True color (ref_mag - mag) of the contaminants.
    prob_tail : float
        Probability that a contaminant passes the cut.
    sigma_c : float
        Intrinsic Gaussian scatter (mag) of the contaminant color. Requires
        z*k < 1, i.e. sigma_c < 2.5*log10(2)/z (~0.20 for prob_tail = 1e-4).
    rescale_factor : float
        Multiplier on the flux errors of all bands.
    include_color_lim : bool
        If True, also return the equivalent color threshold.

    Returns
    -------
    f_ref_lim : array
        Upper limit on f_ref (nanomaggies). Compare against the unclipped
        continuum flux, not one derived from a clipped ref_mag.
    color_lim : array, only if include_color_lim is True
        Equivalent lower limit on color = ref_mag - mag; select
        color > color_lim. It is inf where f_ref_lim <= 0 < f_b (no object
        with positive f_ref passes) and NaN where f_b <= 0 (color undefined).
    '''
    f_b = np.asarray(f_b, dtype=float)
    z = norm.isf(prob_tail)  # gaussian sigma at tail probability prob_tail
    a = 10**(0.4 * c0)
    sigma_b = rescale_factor * 10**(-0.4 * (depth - 22.5)) / 5  # error on selection-band flux
    sigma_x = sigma_b * np.sqrt(1 + a**2 / 4)  # error on f_b - a*f_ref
    # solve f_b - y = z * sqrt(sigma_x**2 + (k*y)**2) for y = a * f_ref_lim (root with y < f_b)
    k = (10**(0.4 * z * sigma_c) - 1) / z
    A = 1 - (z * k)**2
    if A <= 0:
        raise ValueError('sigma_c too large for this prob_tail')
    disc = f_b**2 - A * (f_b**2 - (z * sigma_x)**2)
    f_ref_lim = (f_b - np.sqrt(disc)) / A / a

    if not include_color_lim:
        return f_ref_lim
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            color_lim = np.where(f_ref_lim > 0, 2.5 * np.log10(f_b / f_ref_lim), np.inf)
        color_lim = np.where(f_b > 0, color_lim, np.nan)
        return f_ref_lim, color_lim


def flux_to_mag(flux):
    return 22.5 - 2.5 * np.log10(np.clip(flux, fnodet, None))


def print_density(mask, label, n_prev=None):
    '''Print the density within density_rad_deg and the fraction removed (or added) relative to n_prev.
    Returns the number of objects, to be passed as n_prev to the next call.'''
    n = np.sum(mask & mask_rad)
    msg = f'{n / area:.1f}/deg2\t{label}'
    if n_prev:
        frac = 1 - n / n_prev
        msg += f' ({frac:.2%} removed)' if frac >= 0 else f' ({-frac:.2%} added)'
    print(msg)
    return n


def select_lae(cat, band, cfg):
    print(f'\n{band}')
    mask = np.full(len(cat), True)
    n = print_density(mask, 'initial')

    # cut on maskbits
    for bit in cut_bits:
        mask &= (cat['maskbits'] & 2**bit) == 0
    n = print_density(mask, 'after maskbits cut', n)

    # cut on IBIS fibmag
    fibmag = cat[f'{band}fibmag']
    mask &= (fibmag > FIBMAG_MIN) & (fibmag < FIBMAG_MAX)
    n = print_density(mask, f'after {FIBMAG_MIN} < fibmag < {FIBMAG_MAX} cut', n)

    # cut on excess IBIS flux in the selection band
    # f_ref_lim = flux_threshold(flux_ec[band], depth=median_depths[band], c0=contam_median_colors[band],
    f_ref_lim = flux_threshold(flux_ec[band], depth=cat[band+'_depth_ebv'], c0=contam_median_colors[band],
                                prob_tail=PROB_TAIL, sigma_c=cfg['color_offset_sigma'])
    mask &= cat[f'{band}_ref_flux'] < 10**(-cfg['color_offset'] / 2.5) * f_ref_lim
    n = print_density(mask, 'after IBIS color cut', n)

    # cut on g-r and r-i
    mask1 = (cat['gmag'] - cat['rmag'] < cfg['gr_colmax']) | ((cat['gmag'] < magnodet) & (cat['rmag'] == magnodet)) | ((cat['gmag'] == magnodet) & (cat['rmag'] == magnodet))
    n = print_density(mask & mask1, f'after g-r < {cfg["gr_colmax"]} cut', n)
    mask1 &= (cat['rmag'] - cat['imag'] < cfg['ri_colmax']) | ((cat['rmag'] < magnodet) & (cat['imag'] == magnodet)) | ((cat['rmag'] == magnodet) & (cat['imag'] == magnodet))
    n = print_density(mask & mask1, f'after r-i < {cfg["ri_colmax"]} cut', n)

    # add back rmag>RMAG_ADD_BACK_LIMIT sources that failed the g-r and r-i cuts
    mask1 |= cat['rmag'] > RMAG_ADD_BACK_LIMIT
    mask &= mask1
    n = print_density(mask, f'after adding back rmag>{RMAG_ADD_BACK_LIMIT} sources', n)

    # remove shreds of largish galaxies and brightish stars
    mask &= cat[f'fracflux_{band}'] < 1
    n = print_density(mask, 'after removing shreds', n)

    # better safe than sorry... remove a handful of dubious bright r-mag targets
    bright = cat['rmag'] < RMAG_BRIGHT_CUT
    print(f'remove {np.sum(mask & bright)} targets with r_mag < {RMAG_BRIGHT_CUT}')
    mask &= ~bright

    return mask


##############################################################################################################################
# Load catalogs

cut_bits = get_maskbits_cut_bits(ibis_filters, maskbits_header_fn)

ibis_columns = ['ibis_id_dr1', 'ra', 'dec', 'ebv', 'maskbits']
ibis_columns += [f'{prefix}_{band}' for band in ibis_filters for prefix in ['flux', 'fiberflux', 'fracflux', 'psfdepth']]
hsc_columns = ['ibis_id_dr1'] + [f'forced_flux_{band}' for band in hsc_filters]

ibis = Table(fitsio.read(ibis_fn, columns=ibis_columns))
hsc = Table(fitsio.read(hsc_fn, columns=hsc_columns))
# ibis = Table(fitsio.read(ibis_fn))
# hsc = Table(fitsio.read(hsc_fn))
assert np.all(ibis['ibis_id_dr1'] == hsc['ibis_id_dr1'])

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
print('radius cut', len(cat), np.sum(mask) / len(mask))

mask_rad = cat['CENTER_OFFSET'] < density_rad_deg

##############################################################################################################################
# Compute magnitudes

# extinction-corrected IBIS total fluxes
flux_ec = {band: cat[f'flux_{band}'] * 10**(0.4 * ext_coeffs[band] * cat['ebv']) for band in ibis_filters}

for band in ibis_filters:
    cat[band+'mag'] = flux_to_mag(flux_ec[band])
    cat[band+'fibmag'] = flux_to_mag(cat[f'fiberflux_{band}'] * 10**(0.4 * ext_coeffs[band] * cat['ebv']))

for band in hsc_filters:
    # HSC r2 and i2 mags are still rmag and imag
    cat[band[0]+'mag'] = flux_to_mag(cat[f'forced_flux_{band}'] * 10**(0.4 * ext_coeffs[band[0]] * cat['ebv']))

# leave-one-out average flux of the other 4 IBIS bands
for band in ibis_filters:
    ref_flux = np.zeros(len(cat))
    for band1 in ibis_filters:
        if band1 != band:
            ref_flux += flux_ec[band1] / 4
    cat[band+'_ref_flux'] = ref_flux
    cat[band+'_ref_mag'] = flux_to_mag(ref_flux)

# extinction-corrected 5-sigma PSF depths
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    for band in ibis_filters:
        cat[band+'_depth_ebv'] = -2.5 * (np.log10(5 / np.sqrt(cat[f'psfdepth_{band}'])) - 9) - ext_coeffs[band] * cat['ebv']

median_depths = {}
print('Median depths:')
for band in ibis_filters:
    median_depths[band] = np.nanmedian(cat[band+'_depth_ebv'][mask_rad])
    print(band, f'{median_depths[band]:.2f}')

##############################################################################################################################
# Select LAEs

for band, cfg in LAE_CONFIGS.items():
    cat[f'{band}_sel'] = select_lae(cat, band, cfg)

print('\nAll combined:')
cat['lae_sel'] = False
for band in ibis_filters:
    cat['lae_sel'] |= cat[f'{band}_sel']
    band_density = np.sum(cat[f'{band}_sel'] & mask_rad) / area
    print('{:8.1f}/deg2 after adding {} ({:.1f}/deg2)'.format(np.sum(cat['lae_sel'] & mask_rad) / area, band, band_density))

print('\nFlag previously observed objects')
obs = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/misc/observed_in_t44_t49.fits'))
idx1, idx2, d2d, d_ra, d_dec = match_coord(obs['TARGET_RA'], obs['TARGET_DEC'], cat['ra'], cat['dec'], search_radius=0.5, plot_q=False, keep_all_pairs=True)
idx = np.unique(idx2)
cat['new_obs'] = True
cat['new_obs'][idx2] = False
print('New targets fraction (within {:.2f} radius): {:.2f}%'.format(density_rad_deg, 100*np.sum(cat['new_obs'] & mask_rad & cat['lae_sel'])/np.sum(mask_rad & cat['lae_sel'])))
print('New targets:')
cat['lae_sel'] = False
for band in ibis_filters:
    cat['lae_sel'] |= cat[f'{band}_sel']
    band_density = np.sum(cat[f'{band}_sel'] & mask_rad & cat['new_obs']) / area
    print('{:8.1f}/deg2 after adding {} ({:.1f}/deg2)'.format(np.sum(cat['lae_sel'] & mask_rad & cat['new_obs']) / area, band, band_density))

cat.write(output_all_fn, overwrite=False)
cat[cat['lae_sel']].write(output_fn, overwrite=False)
