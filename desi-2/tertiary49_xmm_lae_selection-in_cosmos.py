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

def get_ext_coeffs():

    return {
        "g": 3.240,
        "r": 2.276,
        "i": 1.633,
        "z": 1.263,
        "y": 1.076,
        "M411": 4.290,
        "M438": 4.099,
        "M464": 3.877,
        "M490": 3.634,
        "M517": 3.389,
    }


# [decam-chatter 18711] Suggested MASKBITS for IBIS
# 2024-12-20
def get_maskbits_cut_bits(bands):

    cut_names = ["NPRIMARY", "BRIGHT", "MEDIUM", "GALAXY"]
    for band in bands:
        cut_names += ["SATUR_{}".format(band), "ALLMASK_{}".format(band)]

    # pick one tractor file, to grab header
    fn = "/pscratch/sd/d/dstn/ibis3/tractor/032/tractor-0329m042.fits"
    hdr = fitsio.read_header(fn, 0)
    keys = np.array([key for key in hdr if key[:5] == "MBIT_"])
    names = np.array([hdr[key] for key in keys])

    cut_keys = keys[np.in1d(names, cut_names)]
    cut_bits = np.array([int(key.split("_")[-1]) for key in cut_keys])

    return cut_bits


# AR offset field by +0.3 in dec to avoid shallow m464 region
field_ra, field_dec, field_rad_deg = 150.1, 2.182 + 0.3, 1.8
area = np.pi * field_rad_deg ** 2

fnodet = 0.001
magnodet = 22.5 - 2.5 * np.log10(fnodet) # 30

fn = '/dvs_ro/cfs/cdirs/desi/users/raichoor/ibis/ibis4-cosmos-forced-hsc-wide.fits'
cat = Table(fitsio.read(fn))

# remove ambiguous columns (as those MAG do not have ebv-corr,
#   and we use different values here; also remove all x-match with desi infos)
black_keys = []
black_keys += [_ for _ in cat.colnames if _[:4] == "MAG_"]
black_keys += [_ for _ in cat.colnames if _[:7] == "FIBMAG_"]
black_keys += [_ for _ in cat.colnames if "DESI" in _]
print("remove these columns: {}".format(",".join(black_keys)))
cat.remove_columns(black_keys)

# clip the magnitudes at 30.0
for band in ['M411', 'M438', 'M464', 'M490', 'M517']:
    cat[band+'fibmag'] = 22.5 - 2.5*np.log10(np.clip(cat["FIBERFLUX_{}".format(band)]*10**(0.4*get_ext_coeffs()[band]*cat['EBV']), fnodet, None))
for band in ['g', 'r', 'i']:
    if band=='g':
        cat[band+'mag'] = 22.5 - 2.5*np.log10(np.clip(cat["FORCED_FLUX_{}".format(band.upper())]*10**(0.4*get_ext_coeffs()[band]*cat['EBV']), fnodet, None))
    else:
        cat[band+'mag'] = 22.5 - 2.5*np.log10(np.clip(cat["FORCED_FLUX_{}2".format(band.upper())]*10**(0.4*get_ext_coeffs()[band]*cat['EBV']), fnodet, None))

# cut on radius
field_c = SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg, frame="icrs")                                                                                                   
cs = SkyCoord(ra=cat["RA"] * u.deg, dec=cat["DEC"] * u.deg, frame="icrs")
cat["FIELD_RADIUS_OFFSET"] = cs.separation(field_c).deg
mask = cat["FIELD_RADIUS_OFFSET"] < field_rad_deg
print("{:.1f}k/deg2\tafter radius cut".format(mask.sum() / 1000. / area))
cat = cat[mask]
print(len(cat), np.sum(mask)/len(mask))

##################### Add HSC-deep photometry where available #####################
hsc = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-dud-rev-no_mag_limit-reduced.fits', columns=['ra', 'dec']))
print(len(hsc))
mask = (hsc['ra']>100) & (hsc['ra']<200)
idx = np.where(mask)[0]
hsc = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-dud-rev-no_mag_limit-reduced.fits', rows=idx))
print(len(hsc))

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))  # https://github.com/rongpu/Python/blob/master/user_modules/match_coord.py
from match_coord import match_coord
idx1, idx2, d2d, d_ra, d_dec = match_coord(hsc['ra'], hsc['dec'], cat['RA'], cat['DEC'], search_radius=1., plot_q=False)

# Replace the HSC-wide forced photometry with HSC-deep cross-matched photometry
for band in ['g', 'r', 'i']:
    mask = ~np.isnan(hsc['{}_cmodel_flux'.format(band)][idx1])  # do not use the HSC-deep value if it is NaN
    cat[band+'mag'][idx2[mask]] = 22.5 - 2.5*np.log10(np.clip(hsc['{}_cmodel_flux'.format(band)][idx1[mask]]/3630.78*10**(0.4*get_ext_coeffs()[band]*cat['EBV'][idx2[mask]]), fnodet, None))
    cat['hsc_deep_'+band] = False
    cat['hsc_deep_'+band][idx2[mask]] = True

##################### LAE selections #####################

######################## M411 ########################

print('\nM411')

fibmagmin, fibmagmax = 22.0, 25.25
mb_colmin = 0.1
ri_colmax = 0.35
gr_colmax = 0.7

mask = np.full(len(cat), True)

# cut on maskbits
cut_bits = get_maskbits_cut_bits(["M411", "M438"])
for bit in cut_bits:
    mask &= (cat["MASKBITS"] & 2**bit) == 0
print("{:.1f}k/deg2\tafter maskbits cut".format(mask.sum() / 1000. / area))

# cut on IBIS fibmag
mask &= (cat['M411fibmag'] > fibmagmin) & (cat['M411fibmag'] < fibmagmax)
print("{:.1f}k/deg2\tafter {} < fibmag < {} cut".format(mask.sum() / 1000. / area, fibmagmin, fibmagmax))

# cut on IBIS color
mask &= cat['M438fibmag'] - cat['M411fibmag'] > mb_colmin
print("{:.1f}k/deg2\tafter IBIS color > {} cut".format(mask.sum() / 1000. / area, mb_colmin))

# cut on g-r and r-i
mask1 = (cat['gmag'] - cat['rmag'] < gr_colmax) | ((cat['gmag'] < magnodet) & (cat['rmag'] == magnodet))
print("{:.1f}k/deg2\tafter g-r < {} cut".format((mask & mask1).sum() / 1000. / area, gr_colmax))
mask1 &= (cat['rmag'] - cat['imag'] < ri_colmax) | ((cat['rmag'] < magnodet) & (cat['imag'] == magnodet))
print("{:.1f}k/deg2\tafter r-i < {} cut".format((mask & mask1).sum() / 1000. / area, ri_colmax))
mask1 |= (cat['gmag']>25.0) | (cat['rmag']>25.0) | (cat['imag']>25.0)  # also keep everything fainter than mag=25
mask &= mask1
print("{:.1f}k/deg2\tafter adding back faint sources".format(mask.sum() / 1000. / area))

# better safe than sorry... remove a handful of dubious bright r-mag targets
mask1 = cat['rmag'] < 18
print("remove {} targets with r_mag < 18".format(np.sum(mask & mask1)))
mask &= (~mask1)

cat['M411_sel'] = mask.copy()

######################## M438 ########################

print('\nM438')

fibmagmin, fibmagmax = 22.0, 25.25
mb_colmin = 0.31
ri_colmax = 0.35
gr_colmax = 0.7

mask = np.full(len(cat), True)

# cut on maskbits
cut_bits = get_maskbits_cut_bits(["M411", "M438"])
for bit in cut_bits:
    mask &= (cat["MASKBITS"] & 2**bit) == 0
print("{:.1f}k/deg2\tafter maskbits cut".format(mask.sum() / 1000. / area))

# cut on IBIS fibmag
mask &= (cat['M438fibmag'] > fibmagmin) & (cat['M438fibmag'] < fibmagmax)
print("{:.1f}k/deg2\tafter {} < fibmag < {} cut".format(mask.sum() / 1000. / area, fibmagmin, fibmagmax))

# cut on IBIS color
mask &= cat['M411fibmag'] - cat['M438fibmag'] > mb_colmin
print("{:.1f}k/deg2\tafter IBIS color > {} cut".format(mask.sum() / 1000. / area, mb_colmin))

# cut on g-r and r-i
mask1 = (cat['gmag'] - cat['rmag'] < gr_colmax) | ((cat['gmag'] < magnodet) & (cat['rmag'] == magnodet))
print("{:.1f}k/deg2\tafter g-r < {} cut".format((mask & mask1).sum() / 1000. / area, gr_colmax))
mask1 &= (cat['rmag'] - cat['imag'] < ri_colmax) | ((cat['rmag'] < magnodet) & (cat['imag'] == magnodet))
print("{:.1f}k/deg2\tafter r-i < {} cut".format((mask & mask1).sum() / 1000. / area, ri_colmax))
mask1 |= (cat['gmag']>25.0) | (cat['rmag']>25.0) | (cat['imag']>25.0)  # also keep everything fainter than mag=25
mask &= mask1
print("{:.1f}k/deg2\tafter adding back faint sources".format(mask.sum() / 1000. / area))

# better safe than sorry... remove a handful of dubious bright r-mag targets
mask1 = cat['rmag'] < 18
print("remove {} targets with r_mag < 18".format(np.sum(mask & mask1)))
mask &= (~mask1)

cat['M438_sel'] = mask.copy()

######################## M464 ########################

print('\nM464')

fibmagmin, fibmagmax = 22.0, 25.25
mb_colmin = 0.3
ri_colmax = 0.35
gr_colmax = 0.7

mask = np.full(len(cat), True)

# cut on maskbits
cut_bits = get_maskbits_cut_bits(["M438", "M464"])
for bit in cut_bits:
    mask &= (cat["MASKBITS"] & 2**bit) == 0
print("{:.1f}k/deg2\tafter maskbits cut".format(mask.sum() / 1000. / area))

# cut on IBIS fibmag
mask &= (cat['M464fibmag'] > fibmagmin) & (cat['M464fibmag'] < fibmagmax)
print("{:.1f}k/deg2\tafter {} < fibmag < {} cut".format(mask.sum() / 1000. / area, fibmagmin, fibmagmax))

# cut on IBIS color
mask &= cat['M438fibmag'] - cat['M464fibmag'] > mb_colmin
print("{:.1f}k/deg2\tafter IBIS color > {} cut".format(mask.sum() / 1000. / area, mb_colmin))

# cut on g-r and r-i
mask1 = (cat['gmag'] - cat['rmag'] < gr_colmax) | ((cat['gmag'] < magnodet) & (cat['rmag'] == magnodet))
print("{:.1f}k/deg2\tafter g-r < {} cut".format((mask & mask1).sum() / 1000. / area, gr_colmax))
mask1 &= (cat['rmag'] - cat['imag'] < ri_colmax) | ((cat['rmag'] < magnodet) & (cat['imag'] == magnodet))
print("{:.1f}k/deg2\tafter r-i < {} cut".format((mask & mask1).sum() / 1000. / area, ri_colmax))
mask1 |= (cat['gmag']>25.0) | (cat['rmag']>25.0) | (cat['imag']>25.0)  # also keep everything fainter than mag=25
mask &= mask1
print("{:.1f}k/deg2\tafter adding back faint sources".format(mask.sum() / 1000. / area))

# better safe than sorry... remove a handful of dubious bright r-mag targets
mask1 = cat['rmag'] < 18
print("remove {} targets with r_mag < 18".format(np.sum(mask & mask1)))
mask &= (~mask1)

cat['M464_sel'] = mask.copy()

######################## M490 ########################

print('\nM490')

fibmagmin, fibmagmax = 22.0, 25.25
mb_colmin = 0.32
ri_colmax = 0.35
gr_colmax = 0.7

mask = np.full(len(cat), True)

# cut on maskbits
cut_bits = get_maskbits_cut_bits(["M464", "M490"])
for bit in cut_bits:
    mask &= (cat["MASKBITS"] & 2**bit) == 0
print("{:.1f}k/deg2\tafter maskbits cut".format(mask.sum() / 1000. / area))

# cut on IBIS fibmag
mask &= (cat['M490fibmag'] > fibmagmin) & (cat['M490fibmag'] < fibmagmax)
print("{:.1f}k/deg2\tafter {} < fibmag < {} cut".format(mask.sum() / 1000. / area, fibmagmin, fibmagmax))

# cut on IBIS color
mask &= cat['M464fibmag'] - cat['M490fibmag'] > mb_colmin
print("{:.1f}k/deg2\tafter IBIS color > {} cut".format(mask.sum() / 1000. / area, mb_colmin))

# cut on g-r and r-i
mask1 = (cat['gmag'] - cat['rmag'] < gr_colmax) | ((cat['gmag'] < magnodet) & (cat['rmag'] == magnodet))
print("{:.1f}k/deg2\tafter g-r < {} cut".format((mask & mask1).sum() / 1000. / area, gr_colmax))
mask1 &= (cat['rmag'] - cat['imag'] < ri_colmax) | ((cat['rmag'] < magnodet) & (cat['imag'] == magnodet))
print("{:.1f}k/deg2\tafter r-i < {} cut".format((mask & mask1).sum() / 1000. / area, ri_colmax))
mask1 |= (cat['gmag']>25.0) | (cat['rmag']>25.0) | (cat['imag']>25.0)  # also keep everything fainter than mag=25
mask &= mask1
print("{:.1f}k/deg2\tafter adding back faint sources".format(mask.sum() / 1000. / area))

# better safe than sorry... remove a handful of dubious bright r-mag targets
mask1 = cat['rmag'] < 18
print("remove {} targets with r_mag < 18".format(np.sum(mask & mask1)))
mask &= (~mask1)

cat['M490_sel'] = mask.copy()

######################## M517 ########################

print('\nM517')

fibmagmin, fibmagmax = 22.0, 25.25
mb_colmin = 0.29
ri_colmax = 0.4
gr_colmax = 0.8

mask = np.full(len(cat), True)

# cut on maskbits
cut_bits = get_maskbits_cut_bits(["M490", "M517"])
for bit in cut_bits:
    mask &= (cat["MASKBITS"] & 2**bit) == 0
print("{:.1f}k/deg2\tafter maskbits cut".format(mask.sum() / 1000. / area))

# cut on IBIS fibmag
mask &= (cat['M517fibmag'] > fibmagmin) & (cat['M517fibmag'] < fibmagmax)
print("{:.1f}k/deg2\tafter {} < fibmag < {} cut".format(mask.sum() / 1000. / area, fibmagmin, fibmagmax))

# cut on IBIS color
mask &= cat['M490fibmag'] - cat['M517fibmag'] > mb_colmin
print("{:.1f}k/deg2\tafter IBIS color > {} cut".format(mask.sum() / 1000. / area, mb_colmin))

# cut on g-r and r-i
mask1 = (cat['gmag'] - cat['rmag'] < gr_colmax) | ((cat['gmag'] < magnodet) & (cat['rmag'] == magnodet))
print("{:.1f}k/deg2\tafter g-r < {} cut".format((mask & mask1).sum() / 1000. / area, gr_colmax))
mask1 &= (cat['rmag'] - cat['imag'] < ri_colmax) | ((cat['rmag'] < magnodet) & (cat['imag'] == magnodet))
print("{:.1f}k/deg2\tafter r-i < {} cut".format((mask & mask1).sum() / 1000. / area, ri_colmax))
mask1 |= (cat['gmag']>25.0) | (cat['rmag']>25.0) | (cat['imag']>25.0)  # also keep everything fainter than mag=25
mask &= mask1
print("{:.1f}k/deg2\tafter adding back faint sources".format(mask.sum() / 1000. / area))

# better safe than sorry... remove a handful of dubious bright r-mag targets
mask1 = cat['rmag'] < 18
print("remove {} targets with r_mag < 18".format(np.sum(mask & mask1)))
mask &= (~mask1)

cat['M517_sel'] = mask.copy()

######################## All combined ########################
print('\nAll combined:')
cat['lae_sel'] = False
for band in ['M411', 'M438', 'M464', 'M490', 'M517']:
    cat['lae_sel'] |= cat['{}_sel'.format(band)]
    print("{:.1f}k/deg2\tafter adding {}".format(np.sum(cat['lae_sel']) / 1000. / area, band))
print()
print(np.sum(cat['lae_sel']), np.sum(cat['lae_sel'])/len(cat))

mask = cat['lae_sel'].copy()
cat = cat[mask]
print(len(cat))
cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/xmm_highz/tertiary49_xmm_lae_targets-in_cosmos.fits', overwrite=False)

