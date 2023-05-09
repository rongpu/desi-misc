from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool
import healpy as hp

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


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
            hp_table[band+'mag_diff_median'][index] = np.nanmedian(dr9[band+'mag_diff'][idx])
            hp_table[band+'mag_diff_mean'][index] = np.nanmean(dr9[band+'mag_diff'][idx])
            hp_table[band+'mag_diff_std'][index] = np.nanstd(dr9[band+'mag_diff'][idx])
            hp_table[band+'mag_diff_nmad'][index] = nmad(dr9[band+'mag_diff'][idx])
            hp_table[band+'mag_n_objects'][index] = np.sum(np.isfinite(dr9[band+'mag_diff'][idx]))

    return hp_table


dr9 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_dr9_south_offsets.fits'))
dr10 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_dr10_offsets.fits'))

idx1, idx2, d2d, d_ra, d_dec = match_coord(dr9['RA'], dr9['DEC'], dr10['RA'], dr10['DEC'], search_radius=0.1, plot_q=False)
dr9 = dr9[idx1]
dr10 = dr10[idx2]

for band in ['g', 'r', 'z']:
    mask = (dr9[band+'_valid']==False) | (dr10[band+'_valid']==False)
    dr9[band+'mag_diff'][mask] = np.nan

n_processes = 64

for nside in [64, 128, 256]:

    npix = hp.nside2npix(nside)

    pix_allobj = hp.pixelfunc.ang2pix(nside, dr9['RA'], dr9['DEC'], lonlat=True)
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
    hp_table.write('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/more/gaia_xp_dr9_south_offset_maps_{}-common_stars.fits'.format(nside), overwrite=True)
