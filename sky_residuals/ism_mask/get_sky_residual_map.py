# Example:
# salloc -N 1 -C cpu -q interactive -t 4:00:00
# python get_sky_residual_map.py 512 south

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp
from multiprocessing import Pool

from desitarget.geomask import sweep_files_touch_hp


# nside = 512
# field = 'south'
# python get_sky_residual_map.py 512 south
nside, field = int(sys.argv[1]), str(sys.argv[2])

nmad = lambda x: 1.4826*np.median(np.abs(x-np.median(x)))

hp_columns = ['nsource_g', 'sky_median_g', 'sky_nmad_g', 'nsource_r', 'sky_median_r', 'sky_nmad_r', 'nsource_z', 'sky_median_z', 'sky_nmad_z']
hp_columns_dtype = ['int32', 'float32', 'float32', 'int32', 'float32', 'float32', 'int32', 'float32', 'float32']

npix = hp.nside2npix(nside)

nside_coarse = 8
npix_coarse = hp.nside2npix(nside_coarse)

sweep_fns = glob.glob('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0/*.fits'.format(field))
sweep_fns = [os.path.basename(ss) for ss in sweep_fns]

healpix_to_sweep_fn = sweep_files_touch_hp(8, None, sweep_fns)[0]
pix_coarse_list = np.arange(npix_coarse)
hp_sweep_counts = np.array([len(healpix_to_sweep_fn[pix_coarse]) for pix_coarse in pix_coarse_list])
mask = hp_sweep_counts>0
pix_coarse_list = pix_coarse_list[mask]
print(len(pix_coarse_list), npix_coarse)

# Shuffle
np.random.seed(12391)
pix_coarse_list = np.random.choice(pix_coarse_list, size=len(pix_coarse_list), replace=False)


def get_sky_residual_in_healpix(pix_coarse):

    sweep_fn_list = healpix_to_sweep_fn[pix_coarse]

    if len(sweep_fn_list)==0:
        return None

    cat = []
    for sweep_fn in sweep_fn_list:

        sweep_path = os.path.join('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0/'.format(field), sweep_fn)

        tmp = Table(fitsio.read(sweep_path, columns=['RA', 'DEC']))
        tmp_pix = hp.ang2pix(nside_coarse, tmp['RA'], tmp['DEC'], nest=True, lonlat=True)
        idx = np.where(tmp_pix==pix_coarse)
        tmp = Table(fitsio.read(sweep_path, rows=idx, columns=['RA', 'DEC', 'MASKBITS', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'NOBS_G', 'NOBS_R', 'NOBS_Z']))

        ext_columns = ['APFLUX_BLOBRESID_G', 'APFLUX_BLOBRESID_R', 'APFLUX_BLOBRESID_Z']
        sweep_ext_fn = sweep_path.replace('/9.0/', '/9.0-extra/').replace('.fits', '-ex.fits')
        tmp_ext = Table(fitsio.read(sweep_ext_fn, rows=idx, columns=ext_columns))
        tmp = hstack([tmp, tmp_ext])

        cat.append(tmp)

    cat = vstack(cat)

    for band in ['g', 'r', 'z']:
        band_cap = band.upper()
        # in units of nanomaggies per sq. arcsec.
        cat['sky_'+band] = (cat['APFLUX_BLOBRESID_'+band_cap][:, -1]-cat['APFLUX_BLOBRESID_'+band_cap][:, -2]) / (np.pi*7**2-np.pi*5**2)

    maskbits = [1, 11, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    mask = mask_clean & ((cat['NOBS_G']>0) | (cat['NOBS_R']>0) | (cat['NOBS_Z']>0))
    cat = cat[mask]
    if len(cat)==0:
        return None

    # cat['pix'] = hp.ang2pix(nside, cat['RA'], cat['DEC'], nest=True, lonlat=True)
    # cat.sort('pix')

    pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], nest=True, lonlat=True)
    pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)
    pixorder = np.argsort(pix_allobj)

    pixcnts = np.insert(pixcnts, 0, 0)
    pixcnts = np.cumsum(pixcnts)

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_unique
    hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_unique, nest=True, lonlat=True)
    for index, hp_column in enumerate(hp_columns):
        hp_table[hp_column] = np.zeros(len(pix_unique), dtype=hp_columns_dtype[index])

    for index in range(len(pix_unique)):
        idx = pixorder[pixcnts[index]:pixcnts[index+1]]
        for band in ['g', 'r', 'z']:
            band_cap = band.upper()
            mask = (cat['NOBS_'+band_cap][idx]>0) & (cat['FLUX_IVAR_'+band_cap][idx]>0)
            if np.sum(mask)>0:
                hp_table['nsource_'+band][index] = np.sum(mask)
                hp_table['sky_median_'+band][index] = np.median(cat['sky_'+band][idx][mask])  # in units of nanomaggies per sq. arcsec.
                hp_table['sky_nmad_'+band][index] = nmad(cat['sky_'+band][idx][mask])  # in units of nanomaggies per sq. arcsec.

    output_path = '/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/tmp/sky_resid_map_{}_{}_{}.fits'.format(nside, field, pix_coarse)

    # if os.path.isfile(output_path):
    #     return None
    # hp_table.write(output_path, overwrite=True)

    return hp_table


n_process = 256
with Pool(processes=n_process) as pool:
    res = pool.map(get_sky_residual_in_healpix, np.arange(npix_coarse), chunksize=1)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

# hp_fns = glob.glob('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/tmp/sky_resid_map_{}_{}_*.fits'.format(nside, field))
# hp_table = []
# for hp_fn in hp_fns:
#     hp_table.append(Table(fitsio.read(hp_fn)))

hp_table = vstack(res)
hp_table['HPXPIXEL'] = hp.nest2ring(nside, hp_table['HPXPIXEL'])
hp_table.sort('HPXPIXEL')
hp_table.write('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_{}_{}.fits'.format(nside, field), overwrite=True)

