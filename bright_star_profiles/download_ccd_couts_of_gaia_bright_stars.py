from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits
import healpy as hp
from astropy import wcs
import tarfile

field = 'south'
region_name = 'decals_ngc'

if (field=='north') and ((band=='g') or (band=='r')):
    pixscale = 0.454
else:
    pixscale = 0.262

download_dir = '/global/homes/r/rongpu/temp/tmp/dr8_ccd_cutouts'
gaia_output_path = '/global/homes/r/rongpu/notebooks/bright_star_profiles/data/gaia_sample_for_ccd_cutouts-decals_ngc.fits'

ramin, ramax, decmin, decmax = 160, 210, -5, 20

######################## Load GAIA catalog ########################

gaia_dir = '/project/projectdirs/cosmo/data/gaia/dr2/healpix'

gaia_nside = 32
gaia_npix = hp.nside2npix(gaia_nside)
gaia_hp_ra, gaia_hp_dec = hp.pix2ang(gaia_nside, np.arange(gaia_npix), nest=True, lonlat=True)
if ramin<ramax:
    gaia_list = np.where((gaia_hp_ra>ramin) & (gaia_hp_ra<ramax) & (gaia_hp_dec>decmin) & (gaia_hp_dec<decmax))[0]
else:
    gaia_list = np.where(((gaia_hp_ra>ramin) | (gaia_hp_ra<ramax)) & (gaia_hp_dec>decmin) & (gaia_hp_dec<decmax))[0]

gaia = []
for hp_idx in gaia_list:
    gaia_fn = (5-len(str(hp_idx)))*'0'+str(hp_idx)
    tmp = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), colnames=['PHOT_G_MEAN_MAG'])
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
    
if ramin<ramax:
    mask = (gaia['ra']>ramin) & (gaia['ra']<ramax) & (gaia['dec']>decmin) & (gaia['dec']<decmax)
else:
    mask = ((gaia['ra']>ramin) | (gaia['ra']<ramax)) & (gaia['dec']>decmin) & (gaia['dec']<decmax)
gaia = gaia[mask]
print(len(gaia))

if field=='north':
    mask = gaia['dec']>32.375
else:
    mask = gaia['dec']<=32.375    
gaia = gaia[mask]
print(len(gaia))
    
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
transform_interp = {}
for bandtmp in ['g', 'r', 'z']:
    if field=='north':
        tmp = Table.read('/global/homes/r/rongpu/notebooks/bright_star_profiles/data/gaia_bassmzls_transform.txt', format='ascii.commented_header')
    else:
        tmp = Table.read('/global/homes/r/rongpu/notebooks/bright_star_profiles/data/gaia_decals_transform.txt', format='ascii.commented_header')
    transform_interp[bandtmp] = interp1d(tmp['bp_rp'], tmp['ls_'+bandtmp], bounds_error=False, fill_value='extrapolate', kind='linear')
    gaia['ls_'+bandtmp] = gaia['phot_g_mean_mag'] + transform_interp[bandtmp](gaia['bp_rp'])

# select 20 stars for every 0.5 magnitude
gaia_g_list = np.arange(6.5, 12.1, 0.5)
print(gaia_g_list)

gaia_idx_all = []
for index in range(len(gaia_g_list)-1):
    gaia_g_min, gaia_g_max = gaia_g_list[index], gaia_g_list[index+1]
    idx = np.where((gaia['phot_g_mean_mag']>gaia_g_min) & (gaia['phot_g_mean_mag']<gaia_g_max))[0]
    if len(idx)>20:
        np.random.seed(123)
        idx = np.random.choice(idx, 20, replace=False)
    gaia_idx_all.append(idx)
    
gaia_idx_all = np.concatenate(gaia_idx_all)
gaia = gaia[gaia_idx_all]

if not os.path.isfile(gaia_output_path):
    gaia.write(gaia_output_path)

# doanload cutouts
for gaia_index in range(len(gaia)):
    
    print(gaia_index, 'phot_g_mean_mag =', gaia['phot_g_mean_mag'][gaia_index])
    ra, dec = gaia['ra'][gaia_index], gaia['dec'][gaia_index]

    size_str = '200' # 200 is the maximum cutout size

    file_path = os.path.join(download_dir, field, '{}_{}.tgz'.format(gaia['source_id'][gaia_index], size_str))
    if not os.path.exists(os.path.dirname(file_path)):
        os.makedirs(os.path.dirname(file_path))

    if (not os.path.isfile(file_path)) or (os.stat(file_path).st_size==0):
        url = 'http://legacysurvey.org/viewer/cutouts-tgz/?ra={:f}&dec={:f}&size={}&layer=dr8-{}'.format(ra, dec, size_str, field)
        cmd = 'wget -O '+file_path+' \"'+url+'\"'
        print(cmd)
        os.system(cmd)

    tar = tarfile.open(file_path)
    img_dir = os.path.join(os.path.dirname(file_path), tar.getnames()[0])
    if not os.path.exists(img_dir):
        tar.extractall(path=os.path.dirname(file_path))