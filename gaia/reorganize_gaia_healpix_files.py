# Convert the files from .csv.gz to .fits

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

import healpy as hp
from multiprocessing import Pool

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/xp_synthetic_decam_photometry/'
phot_fns = sorted(glob.glob('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/xp_synthetic_decam_photometry/hardtouseformat/XpSyntheticDECam_*.fits'))
print(len(phot_fns))

hp_firsts, hp_lasts = [], []
for fn in phot_fns:
    hp_first, hp_last = [int(ii) for ii in os.path.basename(fn).replace('XpSyntheticDECam_', '').replace('.fits', '').split('-')]
    # print(hp_first, hp_last)
    hp_firsts.append(hp_first)
    hp_lasts.append(hp_last)
hp_firsts, hp_lasts = np.array(hp_firsts), np.array(hp_lasts)

gaia_fns = sorted(glob.glob('/global/cfs/cdirs/cosmo/data/gaia/edr3/healpix/healpix-*.fits'))


def do_something(fn):

    gaia = Table(fitsio.read(fn, columns=['SOURCE_ID', 'RA', 'DEC']))

    hp_32 = int(os.path.basename(fn).replace('healpix-', '').replace('.fits', ''))
    npix = hp.nside2npix(32)
    tmp = np.full(npix, False)
    tmp[hp_32] = True
    hp_256_list = np.where(hp.ud_grade(tmp, 256, order_in='NESTED', order_out='NESTED'))[0]

    mask = np.full(len(phot_fns), False)
    for hp_256 in hp_256_list:
        mask |= (hp_256>=hp_firsts) & (hp_256<=hp_lasts)

    file_idx = np.where(mask)[0]
    print(len(file_idx))

    cat = []
    for ii in file_idx:
        cat.append(Table(fitsio.read(phot_fns[ii])))
    cat = vstack(cat)

    mask = np.in1d(cat['source_id'], gaia['SOURCE_ID'])
    cat = cat[mask]

    mask = np.in1d(gaia['SOURCE_ID'], cat['source_id'])
    gaia = gaia[mask]

    if len(gaia)!=len(cat) or not np.all(np.unique(gaia['SOURCE_ID'])==np.unique(cat['source_id'])):
        raise ValueError('gaia and cat have different source_id list')

    if not all(cat['source_id']==gaia['SOURCE_ID']):
        # Matching cat to gaia
        t1_reverse_sort = np.array(gaia['SOURCE_ID']).argsort().argsort()
        cat = cat[np.argsort(cat['source_id'])[t1_reverse_sort]]

    # cat['ra'] = gaia['RA']
    # cat['dec'] = gaia['DEC']

    cat.write(os.path.join(output_dir, 'XpSyntheticDECam_{}.fits'.format(str(hp_32).zfill(5))))
    gaia.write(os.path.join(output_dir, 'Gaia_EDR3/XpSyntheticDECam_{}_Gaia_EDR3.fits'.format(str(hp_32).zfill(5))))

    return len(cat)


n_process = 256
with Pool(processes=n_process) as pool:
    res = pool.map(do_something, gaia_fns, chunksize=1)

counter = np.sum(res)
print('Total number of objects:', counter)

# Sanity check that no object is lost
counter1 = 0
for fn in phot_fns:
    counter1 += fitsio.read_header(fn, ext=1)['NAXIS2']
print(counter==counter1, counter1, counter1-counter)

# Result:
# Total numbert of objects in new files: 219120650
# Total numbert of objects in original files: 219125259
# 4609 objects are missing

