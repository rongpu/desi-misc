# In DECaLS

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
import healpy as hp

gaia_output_path = '/global/cfs/cdirs/desi/users/rongpu/misc/gaia_dr2_for_dr10_i_band_halo.fits'

cat = Table(fitsio.read('/global/cscratch1/sd/rongpu/temp/dr10_i_band_photom-trim.fits', columns=['ra', 'dec']))

nside = 32
pix = np.unique(hp.ang2pix(nside, cat['ra'], cat['dec'], nest=True, lonlat=True))

######################## Load GAIA catalog ########################

gaia_dir = '/project/projectdirs/cosmo/data/gaia/dr2/healpix'

gaia = []
for hp_idx in pix:
    gaia_fn = (5-len(str(hp_idx)))*'0'+str(hp_idx)
    tmp = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), columns=['PHOT_G_MEAN_MAG'])
    # select GAIA_G<13 objects
    idx = np.where(tmp['PHOT_G_MEAN_MAG']<13)[0]
    tmp = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), rows=idx)
    tmp = Table(tmp)
    gaia.append(tmp)
gaia = vstack(gaia)
print(len(gaia))

# convert column names to lower case
for col in gaia.colnames:
    gaia.rename_column(col, col.lower())

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

print('astrometric error cut:')
gaia['pm'] = np.sqrt(gaia['pmra']**2 + gaia['pmdec']**2)
gaia['pmerr'] = np.sqrt(gaia['pmra_error']**2 + gaia['pmdec_error']**2)
mask = gaia['pm']<30
print(np.sum(mask)/len(mask))
mask &= gaia['pmerr']<15
print(np.sum(mask)/len(mask))
mask &= gaia['astrometric_excess_noise']==0
print(np.sum(mask)/len(mask))
# plt.hist(gaia['phot_g_mean_mag'], 100, range=(5, 14), log=True, alpha=0.5)
# plt.hist(gaia['phot_g_mean_mag'][mask], 100, range=(5, 14), log=True, alpha=0.5)
# plt.show()

gaia = gaia[mask]
print(len(gaia))

# plt.figure(figsize=(10, 5))
# plt.plot(gaia['ra'][::2], gaia['dec'][::2], '.', ms=0.5)
# plt.show()

######################## Load GAIA catalog ########################

# Apply GAIA-LS transformation

from scipy.interpolate import interp1d

coeffs = dict(
    i = [0.3396481660, -0.6491867119, -0.3330769819, 0.4381097294,
        0.5752125977, -1.4746570523, 1.2979140762, -0.6371018151,
        0.1948940062, -0.0382055596, 0.0046907449, -0.0003296841,
        0.0000101480],)

bands = ['i']    
for i, b in enumerate(bands):
    mag = np.copy(gaia['phot_g_mean_mag'])
    for order, c in enumerate(coeffs[b]):
        mag += c * (gaia['bp_rp'])**order
    gaia['ls_'+b] = mag

# # select 20 stars for every 0.5 magnitude
# gaia_g_list = np.arange(6.5, 12.1, 0.5)
# print(gaia_g_list)

# gaia_idx_all = []
# for index in range(len(gaia_g_list)-1):
#     gaia_g_min, gaia_g_max = gaia_g_list[index], gaia_g_list[index+1]
#     idx = np.where((gaia['phot_g_mean_mag']>gaia_g_min) & (gaia['phot_g_mean_mag']<gaia_g_max))[0]
#     if len(idx)>20:
#         np.random.seed(123)
#         idx = np.random.choice(idx, 20, replace=False)
#     gaia_idx_all.append(idx)
    
# gaia_idx_all = np.concatenate(gaia_idx_all)
# gaia = gaia[gaia_idx_all]

if not os.path.isfile(gaia_output_path):
    gaia.write(gaia_output_path)