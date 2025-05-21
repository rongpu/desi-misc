# Add stellar density from Yifei's brick list

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr11-early/survey-ccds-decam-dr11-merged-incl-early.fits'
ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'ccdname', 'ra', 'dec']))
print(len(ccd))

bricks = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/bricklist_from_yifei/bricks-exist-ls-dr11-early-gaiagbt19p5.fits', ext=1))

idx1, idx2, d2d, d_ra, d_dec = match_coord(bricks['ra'], bricks['dec'], ccd['ra'], ccd['dec'], search_radius=900., plot_q=False, keep_all_pairs=True)

ccd = ccd[idx2]
print(len(ccd))

ccd['brickid'] = bricks['brickid'][idx1]
ccd['density'] = bricks['density'][idx1]

def get_mean_density(expnum):
    mask = ccd['expnum']==expnum
    return np.mean(ccd['density'][mask])

expnum_all = np.unique(ccd['expnum'])

n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(get_mean_density, expnum_all)
tt = Table()
tt['expnum'] = expnum_all
tt['density_exp_mean'] = np.array(res)
ccd = join(ccd, tt, keys='expnum', join_type='left')
print(len(ccd))

ccd.write('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/survey-ccds-decam-dr11-brickid_stellar_density.fits')
