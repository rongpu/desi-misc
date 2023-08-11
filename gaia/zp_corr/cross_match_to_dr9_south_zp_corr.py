from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

import healpy as hp
from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/dr9_south_cross_match_zp_corr'

sweep_fns = sorted(glob.glob('/pscratch/sd/r/rongpu/dr9_desi_photoz/sweep_zp_corrected/sweep-*.fits'))


def do_something(fn):

    sweep_output_fn = os.path.join(output_dir, os.path.basename(fn).replace('.fits', '-ls.fits'))
    gaia_output_fn = os.path.join(output_dir, os.path.basename(fn).replace('.fits', '-gaia.fits'))

    if os.path.isfile(sweep_output_fn) and os.path.isfile(gaia_output_fn):
        return None

    cat = Table(fitsio.read(fn, columns=['REF_CAT']))
    idx = np.where(cat['REF_CAT']=='G2')[0]
    cat = Table(fitsio.read(fn, rows=idx))

    mask = cat['GAIA_DUPLICATED_SOURCE']==False
    mask &= cat['TYPE']!='DUP'
    mask &= cat['TYPE']=='PSF'
    cat = cat[mask]

    if len(cat)==0:
        return None

    pix_list = np.unique(hp.ang2pix(32, cat['RA'], cat['DEC'], nest=True, lonlat=True))
    gaia = []
    for pix in pix_list:
        tmp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/xp_synthetic_decam_photometry/XpSyntheticDECam_{}.fits'.format(str(pix).zfill(5))))
        tmp1 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/xp_synthetic_decam_photometry/Gaia_EDR3/XpSyntheticDECam_{}_Gaia_EDR3.fits'.format(str(pix).zfill(5))))
        if not all(tmp['source_id']==tmp1['SOURCE_ID']):
            raise ValueError
        tmp.remove_column('source_id')
        tmp = hstack([tmp1, tmp])
        gaia.append(tmp)
    gaia = vstack(gaia)

    mask = gaia['DUPLICATED_SOURCE']==False
    gaia = gaia[mask]
    # print(len(gaia))

    base_fn = os.path.basename(fn)
    ramin = int(base_fn[6:9])
    decmin = int(base_fn[10:13]) if base_fn[9]=='p' else -int(base_fn[10:13])
    gaia_in_brick = (gaia['RA']>ramin) & (gaia['RA']<ramin+10) & (gaia['DEC']>decmin) & (gaia['DEC']<decmin+5)

    # Convert from EDR3's 2016.0 to DR2's 2015.5
    gaia['RA_J2015.5'] = gaia['RA'] - 0.5 * gaia['PMRA'] * 1e-3/3600 / np.cos(np.radians(gaia['DEC']))
    gaia['DEC_J2015.5'] = gaia['DEC'] - 0.5 * gaia['PMDEC'] * 1e-3/3600

    idx1, idx2, d2d, d_ra, d_dec = match_coord(cat['RA'], cat['DEC'], gaia['RA_J2015.5'], gaia['DEC_J2015.5'], search_radius=0.1, plot_q=False)
    print(len(idx2), np.sum(gaia_in_brick), len(idx2)/np.sum(gaia_in_brick))
    # Preserve the order in the sweep catalog
    ii = np.argsort(idx1)
    idx1 = idx1[ii]
    idx2 = idx2[ii]
    cat = cat[idx1]
    gaia = gaia[idx2]

    cat.write(sweep_output_fn, overwrite=True)
    gaia.write(gaia_output_fn, overwrite=True)

    return None


n_process = 128
with Pool(processes=n_process) as pool:
    res = pool.map(do_something, sweep_fns, chunksize=1)


