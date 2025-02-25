# Compute the PSFEx first moments and max location to identify outliers
# srun -N 1 -C cpu -c 256 -t 04:00:00 -q interactive python compute_psfex_stats.py > compute_psfex_stats.log

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

from multiprocessing import Pool
from scipy import ndimage

time_start = time.time()

pixscale = 0.262  # arcsec

surveyccd_path = '/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/dr11/survey-ccds-decam-dr11-merged.fits'
psfex_dir = '/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/dr11/calib/psfex'

# gaia_dir = '/dvs_ro/cfs/cdirs/cosmo/data/gaia/dr3/healpix'


def psfex_stats(ccd_index):
    image_filename = ccd['image_filename'][ccd_index]
    psfex_filename = image_filename.replace('.fits.fz', '-psfex.fits')
    psfex_path = os.path.join(psfex_dir, psfex_filename)
    data = Table(fitsio.read(psfex_path, ext=1))

    mask = data['ccdname']==ccd['ccdname'][ccd_index]
    # if np.sum(mask)!=1:
    #     print(np.sum(mask), 'CCDs match the ccdname')
    #     return 0, 0, 0, 0, 0, 0
    index = np.where(mask)[0][0]
    data = data[index]

    psf = np.array(data['psf_mask'][0])

    # centroid location
    y_cent, x_cent = ndimage.center_of_mass(psf)
    # maximum pixel location
    y_max, x_max = np.unravel_index(np.argmax(psf), psf.shape)

    pix_center = 0.5*(psf.shape[0]-1)
    # flux ratio between center pixel and the pixel with maximum flux
    center_max_ratio = psf[int(np.round(pix_center)), int(np.round(pix_center))] / psf[y_max, x_max]

    # flux ratio between centroid pixel and the pixel with maximum flux
    try:
        centroid_max_ratio = psf[int(np.round(y_cent)), int(np.round(x_cent))] / psf[y_max, x_max]
    except:
        centroid_max_ratio = np.nan
        print('Error: expnum = {}'.format(ccd['expnum'][ccd_index]))
        print(y_cent, x_cent, y_max, x_max)

    # set center pixel location to 0, 0
    x_cent, y_cent, x_max, y_max = x_cent-pix_center, y_cent-pix_center, x_max-pix_center, y_max-pix_center

    ################ Compute ellipticity################
    yy, xx = np.indices(psf.shape).astype('float') - pix_center
    psf1 = psf.copy()
    mask = np.sqrt(xx**2 + yy**2) * pixscale > 4.  # remove the noise-dominated outer part of the model from the ellipticity calculation
    psf1[mask] = 0.
    xx -= x_cent
    yy -= y_cent
    # Compute second moments
    I = psf1 / psf1.sum()  # Normalize intensity
    Qxx = np.sum(xx**2 * I)
    Qyy = np.sum(yy**2 * I)
    Qxy = np.sum(xx * yy * I)
    # Compute ellipticity components
    e1 = (Qxx - Qyy) / (Qxx + Qyy)
    e2 = (2 * Qxy) / (Qxx + Qyy)
    ellipticity = np.sqrt(e1**2 + e2**2)

    # noise equivalent area (in pixels)
    nea = np.sum(psf)**2/np.sum(psf**2)

    fwhm_psfex = data['psf_fwhm'] * pixscale  # in arcsec

    return x_cent, y_cent, x_max, y_max, ellipticity, centroid_max_ratio, center_max_ratio, nea, fwhm_psfex


ccd = Table(fitsio.read(surveyccd_path, columns=['filter', 'ccdname', 'expnum', 'image_filename', 'ccd_cuts', 'image_hdu', 'ra', 'dec', 'propid', 'fwhm']))
print(len(ccd))

# Remove potentially bad CCDs
ccd_blacklist = ['N10', 'N13', 'N15', 'N30', 'S30', 'S7 ']
mask = np.in1d(ccd['ccdname'], ccd_blacklist)
ccd = ccd[~mask]
print(len(ccd))

# Unique exposures
print(len(ccd), len(np.unique(ccd['expnum'])))
# prefer N1, N7, S1:
ccd['priority'] = 100
for ii, ccdname in enumerate(['N1', 'N7', 'S1']):
    mask = ccd['ccdname']==ccdname
    ccd['priority'][mask] = ii
ccd.sort('priority')
_, idx_keep = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx_keep]
ccd.remove_column('priority')
print(len(ccd), len(np.unique(ccd['expnum'])))

n_process = 256
with Pool(processes=n_process) as pool:
    res = pool.map(psfex_stats, np.arange(len(ccd)))
res = np.array(res)

ccd['x_cent'] = res[:, 0]
ccd['y_cent'] = res[:, 1]
ccd['x_max'] = res[:, 2]
ccd['y_max'] = res[:, 3]
ccd['ellipticity'] = res[:, 4]
ccd['centroid_max_ratio'] = res[:, 5]
ccd['center_max_ratio'] = res[:, 6]
ccd['nea'] = res[:, 7]
ccd['fwhm_psfex'] = res[:, 8]

ccd['fwhm_cp'] = ccd['fwhm'] * pixscale  # in arcsec
ccd['fwhm_nea'] = np.sqrt(ccd['nea'] / (4 * np.pi)) * 2.3548 * pixscale  # NEA to FWHM (a la tractor); in arcsec
ccd['r_max'] = np.sqrt((ccd['x_max'])**2+(ccd['y_max'])**2)
ccd['r_cent'] = np.sqrt(ccd['x_cent']**2 + ccd['y_cent']**2)

ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11/survey-ccds-decam-dr11-psfex-stats.fits', overwrite=True)

print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))
