from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool
import healpy as hp

nmad = lambda x: 1.4826 * np.nanmedian(np.abs(x-np.nanmedian(x)))


def get_systematics(pix_idx):

    pix_list = pix_unique[pix_idx]

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)
    for band in ['g', 'r', 'z']:
        hp_table[band+'mag_diff_median'] = np.nan
        hp_table[band+'mag_diff_mean'] = np.nan
        hp_table[band+'mag_diff_std'] = np.nan
        hp_table[band+'mag_diff_nmad'] = np.nan
        hp_table[band+'mag_n_objects'] = 0

    for index in np.arange(len(pix_idx)):

        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]

        for band in ['g', 'r', 'z']:
            hp_table[band+'mag_diff_median'][index] = np.nanmedian(cat[band+'mag_diff'][idx])
            hp_table[band+'mag_diff_mean'][index] = np.nanmean(cat[band+'mag_diff'][idx])
            hp_table[band+'mag_diff_std'][index] = np.nanstd(cat[band+'mag_diff'][idx])
            hp_table[band+'mag_diff_nmad'][index] = nmad(cat[band+'mag_diff'][idx])
            hp_table[band+'mag_n_objects'][index] = np.sum(np.isfinite(cat[band+'mag_diff'][idx]))

    return hp_table


cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_dr9_south_offsets_zp_corr.fits'))

for band in ['g', 'r', 'z']:
    mask = cat[band+'_valid']==False
    cat[band+'mag_diff'][mask] = np.nan

n_processes = 128

# for nside in [64, 128, 256]:
for nside in [128]:

    npix = hp.nside2npix(nside)

    pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
    pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)

    pixcnts = np.insert(pixcnts, 0, 0)
    pixcnts = np.cumsum(pixcnts)

    pixorder = np.argsort(pix_allobj)

    # split among the Cori processors
    pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_processes)

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        res = pool.map(get_systematics, pix_idx_split)

    hp_table = vstack(res)
    hp_table.sort('HPXPIXEL')
    hp_table.write('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_xp_dr9_south_offset_maps_{}_zp_corr.fits'.format(nside), overwrite=True)
