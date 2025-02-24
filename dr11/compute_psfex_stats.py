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

time_start = time.time()


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

    psfex_fwhm = data['psf_fwhm']
    psf = np.array(data['psf_mask'][0])

    grid = np.linspace(-0.5*(psf.shape[0]-1), 0.5*(psf.shape[0]-1), psf.shape[0])
    xx, yy = np.meshgrid(grid, grid)
    psf = psf/np.nansum(psf)
    if np.sum(np.isnan(psf))>0:
        if np.all(np.isnan(psf)):
            print('All Nan: ', psfex_filename)
        else:
            print(np.sum(np.isnan(psf)), 'NaN in: ', psfex_filename)
    
    x_moment = np.sum(xx*psf)
    y_moment = np.sum(yy*psf)
    x_max, y_max = np.unravel_index(np.argmax(psf), psf.shape)
    x_max -= psf.shape[0]
    y_max -= psf.shape[0]

    # noise equivalent area (in pixels)
    nea = np.sum(psf)**2/np.sum(psf**2)

    return x_moment, y_moment, x_max, y_max, nea, psfex_fwhm


ccd = Table(fitsio.read(surveyccd_path, columns=['filter', 'ccdname', 'expnum', 'image_filename', 'ccd_cuts', 'image_hdu', 'ra', 'dec', 'propid', 'fwhm']))
print(len(ccd))

# Remove potentially bad CCDs
ccd_blacklist = ['N10', 'N13', 'N15', 'N30', 'S30', 'S7 ']
mask = np.in1d(ccd['ccdname'], ccd_blacklist)
ccd = ccd[~mask]
print(len(ccd))

# Unique exposures
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
print(len(ccd))

n_process = 256
with Pool(processes=n_process) as pool:
    res = pool.map(psfex_stats, np.arange(len(ccd)))
res = np.array(res)

ccd['x_moment'] = res[:, 0]
ccd['y_moment'] = res[:, 1]
ccd['x_max'] = res[:, 2]
ccd['y_max'] = res[:, 3]
ccd['nea'] = res[:, 4]
ccd['psfex_fwhm'] = res[:, 5]
ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11/survey-ccds-decam-dr11-psfex-stats.fits', overwrite=True)

print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))
