from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from multiprocessing import Pool

n_processes = 64

sweep_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/north/sweep/9.0'

columns = ['TYPE', 'RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 
           'NOBS_G', 'NOBS_R', 'NOBS_Z', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'MASKBITS']

sweep_fn_list = glob.glob(os.path.join(sweep_dir, '*.fits'))

def find_outliers(sweep_fn):

    tmp = fitsio.read(sweep_fn, columns=['FLUX_G', 'FLUX_R', 'FLUX_Z'])
    mask = (tmp['FLUX_G']<=0) | (tmp['FLUX_R']<=0) | (tmp['FLUX_Z']<=0)
    idx = np.where(mask)[0]
    if len(idx)==0:
        print(os.path.basename(sweep_fn), 0)
        return None

    cat = Table(fitsio.read(sweep_fn, columns=columns, rows=idx))
    
    cat['g_outlier'] = False
    cat['r_outlier'] = False
    cat['z_outlier'] = False

    min_nobs = 1
    mask = (cat['NOBS_G']>=min_nobs) & (cat['NOBS_R']>=min_nobs) & (cat['NOBS_Z']>=min_nobs)
    cat = cat[mask]

    maskbits = [1, 5, 6, 7, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))
    cat = cat[mask_clean]
    # print(len(cat))

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cat['gmag'] = 22.5 - 2.5*np.log10(cat['FLUX_G'])
        cat['rmag'] = 22.5 - 2.5*np.log10(cat['FLUX_R'])
        cat['zmag'] = 22.5 - 2.5*np.log10(cat['FLUX_Z'])
        cat['gfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_G'])
        cat['rfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_R'])
        cat['zfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_Z'])

    mask = (cat['rfibermag']<22.5) & (cat['zfibermag']<22.5)
    cat['g_outlier'][mask] = True
    mask = (cat['gfibermag']<23.5) & (cat['zfibermag']<22.0)
    cat['r_outlier'][mask] = True
    mask = (cat['gfibermag']<23.) & (cat['rfibermag']<22.5)
    cat['z_outlier'][mask] = True

    mask = cat['g_outlier'] | cat['r_outlier'] | cat['z_outlier']

    print(os.path.basename(sweep_fn), np.sum(mask), np.sum(cat['g_outlier']), np.sum(cat['r_outlier']), np.sum(cat['z_outlier']))

    if np.sum(mask)>0:
        cat = cat[mask]
        return cat
    else:
        return None


def main():

    print('Starting')

    with Pool(processes=n_processes) as pool:
        res = pool.map(find_outliers, sweep_fn_list)

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    cat = vstack(res)
    cat.write('/global/cfs/cdirs/desi/users/rongpu/dr9/misc/dr9_north_grz_outliers.fits')


    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

