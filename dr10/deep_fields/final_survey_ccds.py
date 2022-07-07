# Apply NCCD cut and restrict to 2.05 radius

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

from multiprocessing import Pool


half_width = 4094/2 * 0.262 / 3600
half_height = 2046/2 * 0.262 / 3600

radec_limits = [[147.8, 152.5, -0.1, 4.5],
                [38.9, 45.1, -3.3, 2.3],
                [32.2, 38.8, -8.7, -2.3],
                [50, 56.9, -31.5, -24.8],
                [4.7, 12.7, -46.3, -40.7]]
radec_limits_new = []
for radec in radec_limits:
    ramin, ramax, decmin, decmax = radec
    radec_limits_new.append([ramin+0.2/np.cos(np.radians(decmin)), ramax-0.2/np.cos(np.radians(decmin)), decmin+0.2, decmax-0.2])
radec_limits = radec_limits_new

field_names = ['COSMOS', 'S1 & S2', 'X1 & X2 & X3 (XMM-LSS)', 'C1 & C2 & C3', 'E1 & E2']

nccd_threshold = 120
nbricks_threshold = 100
n_processes = 16


bricks1 = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields-cosmos.fits'))
bricks2 = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields-des-sn.fits'))
bricks = vstack([bricks1, bricks2])
print('bricks', len(bricks))

surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1.fits'
ccd = Table(fitsio.read(surveyccd_path))
print('ccd', len(ccd))

bricks_all = bricks.copy()
ccd_all = ccd.copy()

band = 'z'
field_index = 3

step_size = 10


ccd = ccd_all.copy()
mask = ccd['filter']==band
ccd = ccd[mask]
radec = radec_limits[field_index]
ramin, ramax, decmin, decmax = radec
mask = (ccd['ra']>ramin) & (ccd['ra']<ramax) & (ccd['dec']>decmin) & (ccd['dec']<decmax)
ccd = ccd[mask]
print('ccd', len(ccd))

bricks = bricks_all.copy()
mask = (bricks['RA']>ramin) & (bricks['RA']<ramax) & (bricks['DEC']>decmin) & (bricks['DEC']<decmax)
bricks = bricks[mask]
print('bricks', len(ccd))

ccds_to_keep = np.arange(len(ccd))


def get_nccd(ccd_index):
    mask = (bricks['RA']>ccd['ra'][ccd_index]-half_width/np.cos(np.radians(bricks['DEC']))) & (bricks['RA']<ccd['ra'][ccd_index]+half_width/np.cos(np.radians(bricks['DEC']))) \
        & (bricks['DEC']>ccd['dec'][ccd_index]-half_height) & (bricks['DEC']<ccd['dec'][ccd_index]+half_height)
    nccd = np.array(mask, dtype=int)
    return nccd


def get_brickscore(ccd_index):
    mask = (bricks['RA']>ccd['ra'][ccd_index]-half_width/np.cos(np.radians(bricks['DEC']))) & (bricks['RA']<ccd['ra'][ccd_index]+half_width/np.cos(np.radians(bricks['DEC']))) \
        & (bricks['DEC']>ccd['dec'][ccd_index]-half_height) & (bricks['DEC']<ccd['dec'][ccd_index]+half_height)
    brickscore = np.sum(bricks['nccd'][mask]-nccd_threshold)
    return brickscore


with Pool(processes=n_processes) as pool:
    res = pool.map(get_nccd, ccds_to_keep)
nccd = np.sum(res, axis=0)
bricks['nccd'] = nccd
mask = bricks['nccd']>30  # Remove the very shallow bricks
bricks = bricks[mask]

counter = 0
while np.sum(bricks['nccd']>nccd_threshold)>nbricks_threshold:
    counter += 1
    print('Iteration {}'.format(counter))

    with Pool(processes=n_processes) as pool:
        brickscore = np.array(pool.map(get_brickscore, ccds_to_keep))

    ccds_to_delete = np.argsort(brickscore/ccd['efftime'][ccds_to_keep])[-step_size:]
    ccds_to_keep = np.delete(ccds_to_keep, ccds_to_delete)
    print('nccd', len(ccds_to_keep))

    with Pool(processes=n_processes) as pool:
        res = pool.map(get_nccd, ccds_to_keep)
    nccd = np.sum(res, axis=0)
    bricks['nccd'] = nccd
    print('Above-threshold bricks:', np.sum(bricks['nccd']>nccd_threshold))


print(len(ccds_to_keep))

ccd = ccd[ccds_to_keep]
ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/test/tmp/survey-ccds-{}-{}.fits'.format(field_index, band), overwrite=True)

