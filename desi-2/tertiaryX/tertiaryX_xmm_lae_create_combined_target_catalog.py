# To-do: add filler targets

import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

hip = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-ibis_deep.fits'))
n_all = len(hip)
hip = hip[hip['new_obs']]
hip['priority'] = 'high'
print('ibis_deep (high):      {} selected, {} new_obs'.format(n_all, len(hip)))

mip = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-ibis_wide.fits'))
n_all = len(mip)
mip = mip[mip['new_obs']]
mip['priority'] = 'mid'
print('ibis_wide (mid):       {} selected, {} new_obs'.format(n_all, len(mip)))

lop = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-2x_dens.fits'))
n_all = len(lop)
lop = lop[lop['new_obs']]
lop['priority'] = 'low'
print('2x_dens (low):         {} selected, {} new_obs'.format(n_all, len(lop)))

vlo = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiaryX/tertiaryX_xmm_lae_targets-filler.fits'))
n_all = len(vlo)
vlo = vlo[vlo['new_obs']]
vlo['priority'] = 'very low'
print('filler (very low):     {} selected, {} new_obs'.format(n_all, len(vlo)))

hip = hip[['ibis_id_dr1', 'ra', 'dec', 'priority']]
mip = mip[['ibis_id_dr1', 'ra', 'dec', 'priority']]
lop = lop[['ibis_id_dr1', 'ra', 'dec', 'priority']]
vlo = vlo[['ibis_id_dr1', 'ra', 'dec', 'priority']]
cat = vstack([hip, mip, lop, vlo])
n_stacked = len(cat)
print('Stacked: {} new_obs targets in total'.format(n_stacked))

_, idx = np.unique(cat['ibis_id_dr1'], return_index=True)
cat = cat[idx]
print('After removing duplicates (keeping highest priority): {} unique targets ({} duplicates removed)'.format(len(cat), n_stacked - len(cat)))
for p in ['high', 'mid', 'low', 'very low']:
    print('    {:9s}: {}'.format(p, np.sum(cat['priority'] == p)))

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

field_ra, field_dec, field_rad_deg = 35.75, -4.75, 1.5

from astropy.coordinates import SkyCoord
from astropy import units as u

# cut on radius
field_c = SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg, frame='icrs')
cs = SkyCoord(ra=cat['ra'] * u.deg, dec=cat['dec'] * u.deg, frame='icrs')
cat['CENTER_OFFSET'] = cs.separation(field_c).deg
mask = cat['CENTER_OFFSET'] < field_rad_deg
area = np.pi * field_rad_deg ** 2
print('Density = {:.1f} per deg2 within {:.2f} degree radius'.format(np.sum(mask)/area, field_rad_deg))
