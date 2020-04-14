from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits
import healpy as hp

from astropy import units as u
from astropy.coordinates import SkyCoord

field = 'south'
galactic_b_limit = 20.

gaia_output_path = '/global/homes/r/rongpu/data/brightest_gaia_stars-{}.fits'.format(field)

# ######################## Load GAIA catalog ########################

gaia_dir = '/project/projectdirs/cosmo/data/gaia/dr2/healpix'

gaia_nside = 32
gaia_npix = hp.nside2npix(gaia_nside)
gaia_hp_ra, gaia_hp_dec = hp.pix2ang(gaia_nside, np.arange(gaia_npix), nest=True, lonlat=True)
c = SkyCoord(gaia_hp_ra, gaia_hp_dec, unit='deg')
tmp = c.galactic
gaia_hp_l, gaia_hp_b = tmp.l.to_value('deg'), tmp.b.to_value('deg')

# ramin, ramax, decmin, decmax = 130, 265, 35, 75
# if ramin<ramax:
#     gaia_list = np.where((gaia_hp_ra>ramin) & (gaia_hp_ra<ramax) & (gaia_hp_dec>decmin) & (gaia_hp_dec<decmax))[0]
# else:
#     gaia_list = np.where(((gaia_hp_ra>ramin) | (gaia_hp_ra<ramax)) & (gaia_hp_dec>decmin) & (gaia_hp_dec<decmax))[0]

gaia_list = np.where(np.abs(gaia_hp_b)>galactic_b_limit)[0]
print(len(gaia_list))

gaia = []
for hp_idx in gaia_list:
    gaia_fn = (5-len(str(hp_idx)))*'0'+str(hp_idx)
    tmp = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), colnames=['PHOT_G_MEAN_MAG'])
    # select GAIA_G<13 objects
    idx = np.where(tmp['PHOT_G_MEAN_MAG']<10)[0]
    tmp = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), rows=idx)
    tmp = Table(tmp)
    gaia.append(tmp)
gaia = vstack(gaia)
print(len(gaia))

# convert column names to lower case
for col in gaia.colnames:
    gaia.rename_column(col, col.lower())
    
# if ramin<ramax:
#     mask = (gaia['ra']>ramin) & (gaia['ra']<ramax) & (gaia['dec']>decmin) & (gaia['dec']<decmax)
# else:
#     mask = ((gaia['ra']>ramin) | (gaia['ra']<ramax)) & (gaia['dec']>decmin) & (gaia['dec']<decmax)
# gaia = gaia[mask]
# print(len(gaia))

# if field=='north':
#     mask = gaia['dec']>32.375
# else:
#     mask = gaia['dec']<=32.375    
# gaia = gaia[mask]
# print(len(gaia))
    
# Remove the globular clusters and bad regions
mask = (gaia['ra']>198.0) & (gaia['ra']<198.5) & (gaia['dec']>18.0) & (gaia['dec']<18.4)
mask |=  (gaia['ra']>199.0) & (gaia['ra']<199.2) & (gaia['dec']>17.6) & (gaia['dec']<17.8)
mask |=  (gaia['ra']>182.4) & (gaia['ra']<182.65) & (gaia['dec']>18.44) & (gaia['dec']<18.64)
mask |= (gaia['ra']>13.0) & (gaia['ra']<13.4) & (gaia['dec']>-26.7) & (gaia['dec']<-26.45)
# mask |= (gaia['ra']>14.5) & (gaia['ra']<15.5) & (gaia['dec']>-34.0) & (gaia['dec']<-33.2)
gaia = gaia[~mask]
print(len(gaia))

gaia['bp_rp'] = gaia['phot_bp_mean_mag'] - gaia['phot_rp_mean_mag']
gaia['pm'] = np.sqrt(gaia['pmra']**2 + gaia['pmdec']**2)

# Remove duplicates
if len(np.unique(gaia['source_id']))<len(gaia):
    print('Duplicates exist!')
    gaia.sort('source_id')
    mask = gaia['source_id'][1:]==gaia['source_id'][:-1]
    mask = np.concatenate([[False], mask])
    gaia = gaia[~mask]
    
# Remove duplicated_source==True
mask = gaia['duplicated_source'].copy()
if np.sum(mask)>0:
    print('{} objects with duplicated_source==True'.format(np.sum(mask)))
    gaia = gaia[~mask]
    print(len(gaia))

# Remove objects with invalid GAIA color
mask = (~np.isfinite(gaia['bp_rp'])) | (~np.isfinite(gaia['phot_g_mean_mag']))
mask |= (gaia['phot_bp_mean_mag']==0) | (gaia['phot_rp_mean_mag']==0)
if np.sum(mask)>0:
    print('{} objects with invalid bp_rp'.format(np.sum(mask)))
    gaia = gaia[~mask]
    print(len(gaia))

# print('astrometric error cut:')
# gaia['pm'] = np.sqrt(gaia['pmra']**2 + gaia['pmdec']**2)
# gaia['pmerr'] = np.sqrt(gaia['pmra_error']**2 + gaia['pmdec_error']**2)
# mask = gaia['pm']<30
# print(np.sum(mask)/len(mask))
# mask &= gaia['pmerr']<15
# print(np.sum(mask)/len(mask))
# mask &= gaia['astrometric_excess_noise']==0
# print(np.sum(mask)/len(mask))
# # plt.hist(gaia['phot_g_mean_mag'], 100, range=(5, 14), log=True, alpha=0.5)
# # plt.hist(gaia['phot_g_mean_mag'][mask], 100, range=(5, 14), log=True, alpha=0.5)
# # plt.show()
# gaia = gaia[mask]
# print(len(gaia))

# plt.figure(figsize=(10, 5))
# plt.plot(gaia['ra'][::2], gaia['dec'][::2], '.', ms=0.5)
# plt.show()

# Apply GAIA-LS transformation
from scipy.interpolate import interp1d
transform_interp = {}
for bandtmp in ['g', 'r', 'z']:
    if field=='north':
        tmp = Table.read('/global/homes/r/rongpu/notebooks/bright_star_profiles/data/gaia_bassmzls_transform.txt', format='ascii.commented_header')
    else:
        tmp = Table.read('/global/homes/r/rongpu/notebooks/bright_star_profiles/data/gaia_decals_transform.txt', format='ascii.commented_header')
    transform_interp[bandtmp] = interp1d(tmp['bp_rp'], tmp['ls_'+bandtmp], bounds_error=False, fill_value='extrapolate', kind='linear')
    gaia['ls_'+bandtmp] = gaia['phot_g_mean_mag'] + transform_interp[bandtmp](gaia['bp_rp'])

gaia.write(gaia_output_path)
