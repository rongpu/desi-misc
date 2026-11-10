# To-do: add filler targets

import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp
from astropy.coordinates import SkyCoord
from astropy import units as u

field_ra, field_dec, field_rad_deg = 35.75, -4.75, 1.5
field_c = SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg, frame='icrs')
area = np.pi * field_rad_deg ** 2

def density(tab):
    '''Number density (per deg2) of the targets within the field radius.'''
    cs = SkyCoord(ra=tab['ra'] * u.deg, dec=tab['dec'] * u.deg, frame='icrs')
    return np.sum(cs.separation(field_c).deg < field_rad_deg) / area

print('Densities are within {:.2f} degree radius of ({}, {})'.format(field_rad_deg, field_ra, field_dec))

hip = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-ibis_deep.fits'))
n_all, d_all = len(hip), density(hip)
hip = hip[hip['new_obs']]
hip['priority'] = 'high'
print('ibis_deep (high):      {} selected ({:.1f}/deg2), {} new_obs ({:.1f}/deg2)'.format(n_all, d_all, len(hip), density(hip)))

mip = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-ibis_wide.fits'))
n_all, d_all = len(mip), density(mip)
mip = mip[mip['new_obs']]
mip['priority'] = 'mid'
print('ibis_wide (mid):       {} selected ({:.1f}/deg2), {} new_obs ({:.1f}/deg2)'.format(n_all, d_all, len(mip), density(mip)))

lop = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-2x_dens.fits'))
n_all, d_all = len(lop), density(lop)
lop = lop[lop['new_obs']]
lop['priority'] = 'low'
print('2x_dens (low):         {} selected ({:.1f}/deg2), {} new_obs ({:.1f}/deg2)'.format(n_all, d_all, len(lop), density(lop)))

vlo = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-filler.fits'))
n_all, d_all = len(vlo), density(vlo)
vlo = vlo[vlo['new_obs']]
vlo['priority'] = 'very low'
print('filler (very low):     {} selected ({:.1f}/deg2), {} new_obs ({:.1f}/deg2)'.format(n_all, d_all, len(vlo), density(vlo)))

hip = hip[['ibis_id_dr1', 'ra', 'dec', 'priority']]
mip = mip[['ibis_id_dr1', 'ra', 'dec', 'priority']]
lop = lop[['ibis_id_dr1', 'ra', 'dec', 'priority']]
vlo = vlo[['ibis_id_dr1', 'ra', 'dec', 'priority']]
cat = vstack([hip, mip, lop, vlo])
n_stacked = len(cat)
print('Stacked: {} new_obs targets in total ({:.1f}/deg2)'.format(n_stacked, density(cat)))

_, idx = np.unique(cat['ibis_id_dr1'], return_index=True)
cat = cat[idx]
print('After removing duplicates (keeping highest priority): {} unique targets ({} duplicates removed), {:.1f}/deg2'.format(len(cat), n_stacked - len(cat), density(cat)))
for p in ['high', 'mid', 'low', 'very low']:
    print('    {:9s}: {} ({:.1f}/deg2)'.format(p, np.sum(cat['priority'] == p), density(cat[cat['priority'] == p])))

mask = np.isin(cat['ibis_id_dr1'], hip['ibis_id_dr1'])
cat['lae_ibis_deep'] = False
cat['lae_ibis_deep'][mask] = True

mask = np.isin(cat['ibis_id_dr1'], mip['ibis_id_dr1'])
cat['lae_ibis_wide'] = False
cat['lae_ibis_wide'][mask] = True

mask = np.isin(cat['ibis_id_dr1'], lop['ibis_id_dr1'])
cat['lae_ibis_2x_dens'] = False
cat['lae_ibis_2x_dens'][mask] = True

mask = np.isin(cat['ibis_id_dr1'], vlo['ibis_id_dr1'])
cat['lae_ibis_filler'] = False
cat['lae_ibis_filler'][mask] = True

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-combined.fits', overwrite=False)
