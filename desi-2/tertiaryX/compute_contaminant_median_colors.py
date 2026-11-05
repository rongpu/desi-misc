# Compute the median IBIS colors of mag_ref - mag for the bulk of the sample before the IBIS color cut
# Most of these objects are contaminants
# Use XMM-LSS because its IBIS imaging is deeper than in COSMOS

import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

from astropy.coordinates import SkyCoord
from astropy import units as u

params = {'legend.fontsize': 'medium',
          'axes.labelsize': 'medium',
          'axes.titlesize': 'medium',
          'xtick.labelsize': 'medium',
          'ytick.labelsize': 'medium',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']

ext_coeffs = {
'g': 3.240, 'r': 2.276, 'i': 1.633, 'z': 1.263, 'y': 1.076,
'M411': 3.774, 'M438': 3.576, 'M464': 3.366, 'M490': 3.141, 'M517': 2.922}

def get_maskbits_cut_bits(bands):

    cut_names = ['NPRIMARY', 'BRIGHT', 'MEDIUM', 'GALAXY']
    cut_names += [f'{prefix}_{band}' for band in bands for prefix in ['SATUR', 'ALLMASK']]

    # pick one tractor file, to grab header (MBIT_<bit> = '<name>')
    fn = '/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/tractor/032/tractor-0329m042.fits'
    hdr = fitsio.read_header(fn, 0)
    name_to_bit = {hdr[key]: int(key[len('MBIT_'):]) for key in hdr if key.startswith('MBIT_')}

    missing = [name for name in cut_names if name not in name_to_bit]
    if missing:
        raise ValueError(f'maskbit names not found in {fn}: {missing}')

    return np.array(sorted(name_to_bit[name] for name in cut_names))

cut_bits = get_maskbits_cut_bits(ibis_filters)
print(cut_bits)

field_ra, field_dec, field_rad_deg = 35.75, -4.75, 1.8

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

mask = cat['CENTER_OFFSET'] < 1.5
cat = cat[mask]
print('radius cut', len(cat), np.sum(mask)/len(mask))

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
    cat[band+'_ref_flux'] = 0.
    for band1 in ibis_filters:
        if band1!=band:
            cat[band+'_ref_flux'] += cat['flux_'+band1] * 10**(0.4*ext_coeffs[band1]*cat['ebv']) / 4
    cat[band+'_ref_mag'] = 22.5 - 2.5*np.log10(np.clip(cat[band+'_ref_flux'], fnodet, None))

FIBMAG_MIN, FIBMAG_MAX = 22.0, 25.0
RMAG_BRIGHT_CUT = 18
MAGNODET = 30.0
MAGNODET_HSC = 25.0
RI_COLMAX = 0.35
GR_COLMAX = 0.7

contam_median_colors = {}

for band in ibis_filters:
    print(f'\n{band}')
    mask = np.full(len(cat), True)
    print(f'{mask.sum()}\tinitial')

    # cut on maskbits
    for bit in cut_bits:
        mask &= (cat['maskbits'] & 2**bit) == 0
    print(f'{mask.sum()}\tafter maskbits cut')

    # cut on IBIS fibmag
    fibmag = cat[f'{band}fibmag']
    mask &= (fibmag > FIBMAG_MIN) & (fibmag < FIBMAG_MAX)
    print(f'{mask.sum()}\tafter {FIBMAG_MIN} < fibmag < {FIBMAG_MAX} cut')

    # cut on g-r and r-i
    mask1 = (cat['gmag'] - cat['rmag'] < GR_COLMAX) | ((cat['gmag'] < MAGNODET) & (cat['rmag'] == MAGNODET)) | ((cat['gmag'] == MAGNODET) & (cat['rmag'] == MAGNODET))
    print(f'{(mask & mask1).sum()}\tafter g-r < {GR_COLMAX} cut')
    mask1 &= (cat['rmag'] - cat['imag'] < RI_COLMAX) | ((cat['rmag'] < MAGNODET) & (cat['imag'] == MAGNODET)) | ((cat['rmag'] == MAGNODET) & (cat['imag'] == MAGNODET))
    print(f'{(mask & mask1).sum()}\tafter r-i < {RI_COLMAX} cut')
    mask1 |= (cat['rmag'] > MAGNODET_HSC)
    print(f'{(mask & mask1).sum()}\tafter adding back rmag>{MAGNODET_HSC} sources')
    mask &= mask1

    # remove shreds of largish galaxies and brightish stars
    mask &= cat['fracflux_'+band]<1
    print(f'{mask.sum()}\tafter removing shreds')
    
    mask1 = mask.copy() & (cat[f'{band}mag']>24.4) & (cat[f'{band}mag']<24.5)
    median_color = np.nanmedian((cat[band+'_ref_mag'] - cat[band+'mag'])[mask1])
    print('median IBIS color = {:.2f}'.format(median_color))
    
    contam_median_colors[band] = median_color

print()
print({k: round(float(v), 2) for k, v in contam_median_colors.items()})
