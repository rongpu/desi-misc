from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool
import healpy as hp


def get_systematics(pix_list):

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)

    return hp_table


def read_file(fn):
    cat = Table(fitsio.read(fn, columns=['RA', 'DEC', 'DUPLICATED_SOURCE']))
    return cat


fns = sorted(glob.glob('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/xp_synthetic_decam_photometry/Gaia_EDR3/*.fits'))

n_process = 128
with Pool(processes=n_process) as pool:
    res = pool.map(read_file, fns, chunksize=1)
cat = vstack(res)
print(len(cat))

mask = cat["DUPLICATED_SOURCE"]==False
cat = cat[mask]
print(len(cat))

nside = 256
npix = hp.nside2npix(nside)
pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
pix_unique, pix_count = np.unique(pix_allobj, return_counts=True)
hp_table = get_systematics(pix_unique)
hp_table['N_GAIA'] = pix_count
hp_table.write('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_xp_density_map_256.fits', overwrite=True)
