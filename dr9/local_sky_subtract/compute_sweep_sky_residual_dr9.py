# THIS CODE HAS NOT BEEN IMPLEMENTED

from __future__ import division, print_function
import sys, os, glob, gc, warnings, time
import numpy as np
from astropy.io import fits
from astropy.table import Table, vstack, hstack
import fitsio
from multiprocessing import Pool

n_processes = 32
field = 'south'
# field = 'north'

sweep_dir = '/global/cscratch1/sd/landriau/dr9m-partial-sweep/{}/sweep'.format(field)
resid_dir = '/global/cscratch1/sd/rongpu/desi/dr9m-oct26-2020_sweep_apflux_resid/{}'.format(field)

sweep_path_all = np.sort(glob.glob(os.path.join(sweep_dir, '*.fits')))
sweep_fn_all = [os.path.basename(sweep_path_all[ii]) for ii in range(len(sweep_path_all))]

columns = ['RA', 'DEC', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']


def compute_sky_residual(sweep_index):
    
    np.random.seed(123)

    sweep_fn = sweep_fn_all[sweep_index]
    print(sweep_fn)

    fits_tmp = fits.open(os.path.join(sweep_dir, sweep_fn))
    fits_length = fits_tmp[1].header['NAXIS2']
    if fits_length<50:
        return None

    # downsampling
    idx = np.sort(np.random.choice(fits_length, size=fits_length//50, replace=False))

    tmp = Table(fitsio.read(os.path.join(sweep_dir, sweep_fn), columns=columns, rows=idx))
    tmp1 = Table(fitsio.read(os.path.join(resid_dir, sweep_fn[:-5]+'-resid.fits'), rows=idx))
    tmp = hstack([tmp, tmp1])
    
    # Requiring 2+ exposures in each band
    mask = (tmp['NOBS_G']>=2) & (tmp['NOBS_R']>=2) & (tmp['NOBS_Z']>=2) & (tmp['apflux_resid_g'][:, 0]!=-99)
    tmp = tmp[mask]

    if len(tmp)==0:
        return None

    for band in ['g', 'r', 'z']:
        tmp[band+'_sky'] = (tmp['apflux_resid_'+band][:, -1]-tmp['apflux_resid_'+band][:, -2]) / (np.pi*7**2-np.pi*5**2)
    for band in ['w1', 'w2']:
        tmp[band+'_sky_5_7'] = (tmp['apflux_resid_'+band][:, 2]-tmp['apflux_resid_'+band][:, 1]) / (np.pi*7**2-np.pi*5**2)
        tmp[band+'_sky_7_9'] = (tmp['apflux_resid_'+band][:, 3]-tmp['apflux_resid_'+band][:, 2]) / (np.pi*9**2-np.pi*7**2)
        tmp[band+'_sky_9_11'] = (tmp['apflux_resid_'+band][:, 4]-tmp['apflux_resid_'+band][:, 3]) / (np.pi*11**2-np.pi*9**2)

    # blob residual
    for band in ['g', 'r', 'z']:
        tmp[band+'_blobsky'] = (tmp['apflux_blobresid_'+band][:, -1]-tmp['apflux_blobresid_'+band][:, -2]) / (np.pi*7**2-np.pi*5**2)

    # tmp = tmp[['RA', 'DEC', 'MASKBITS', 'g_sky', 'r_sky', 'z_sky', 'g_blobsky', 'r_blobsky', 'z_blobsky']]
    gc.collect()

    return tmp


if __name__ == '__main__':

    print('Start!')

    time_start = time.time()

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(compute_sky_residual, range(len(sweep_fn_all)))

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    stacked = vstack(res)
    print(len(stacked))

    stacked.write('/global/cscratch1/sd/rongpu/temp/dr9_sky_residual_'+field+'.fits', overwrite=True)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
