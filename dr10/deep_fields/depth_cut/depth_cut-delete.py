# Despite the name, it's putting a limit on the number of CCDs (and exposures) rather than depth
# parallel --jobs 7 depth_cut_tasks.txt ; exit
# or
# salloc -N 2 -C cpu -q interactive -t 04:00:00
# srun --wait=0 --ntasks-per-node 1 depth_cut_payload.sh depth_cut_tasks.txt

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import argparse
from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


half_width = 4094/2 * 0.262 / 3600
half_height = 2046/2 * 0.262 / 3600

field_names = ['COSMOS', 'S1 & S2', 'X1 & X2 & X3 (XMM-LSS)', 'C1 & C2 & C3', 'E1 & E2']
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

des_nominal_efftimes = {'g': 239.6,
                        'r': 139.9,
                        'i': 53.4,
                        'z': 19.0,
                        'Y': 4.2}


parser = argparse.ArgumentParser()
parser.add_argument('field_index')
parser.add_argument('band')
args = parser.parse_args()
field_index = int(args.field_index)
band = args.band

# field_index = 0
# band = 'g'


radec = radec_limits[field_index]
ramin, ramax, decmin, decmax = radec

nccd_threshold = 120
nbricks_threshold = 100
step_size = 10

n_processes = 16

bricks1 = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields-cosmos.fits'))
bricks2 = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields-des-sn.fits'))
bricks = vstack([bricks1, bricks2])
print('bricks', len(bricks))

mask = (bricks['RA']>ramin) & (bricks['RA']<ramax) & (bricks['DEC']>decmin) & (bricks['DEC']<decmax)
bricks = bricks[mask]
print('bricks', len(bricks))

surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-defringed.fits'
ccd = Table(fitsio.read(surveyccd_path))
print('ccd', len(ccd))

mask = ccd['filter']==band
ccd = ccd[mask]
print('ccd', len(ccd))

mask = (ccd['ra']>ramin) & (ccd['ra']<ramax) & (ccd['dec']>decmin) & (ccd['dec']<decmax)
ccd = ccd[mask]
print('ccd', len(ccd))

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx]
exp['efftime'] = 10**(0.4*exp['zpt']-9) * exp['exptime'] / (exp['median_ccdskycounts'] * exp['median_psf_fwhm']**2)


############################## Reduce the number of little- and no-dither exposures ##############################

pointing_ra = np.array([150.1166, 54.2743, 54.2743, 52.6484, 34.4757, 35.6645, 36.4500, 42.8200, 41.1944, 7.8744, 9.5000])
pointing_dec = np.array([2.2058, -27.1116, -29.0884, -28.1000, -4.9295, -6.4121, -4.6000, 0.0000, -0.9884, -43.0096, -43.9980])
pointing_names = np.array(['COSMOS', 'C1', 'C2', 'C3', 'X1', 'X2', 'X3', 'S1', 'S2', 'E1', 'E2'])

expnum_list_to_discard = []
exp['weight'] = 1.

# for band in ['g', 'r', 'i', 'z', 'Y']:
for index in range(len(pointing_ra)):

    ra, dec, pointing_name = pointing_ra[index], pointing_dec[index], pointing_names[index]

    # Give half the weight for <1 arcsec exposures
    search_radius = 1.
    idx1, idx2, d2d, d_ra, d_dec = match_coord([ra], [dec], exp['ra_bore'], exp['dec_bore'], search_radius=search_radius, plot_q=False, keep_all_pairs=True, verbose=False)
    exp['weight'][idx2] = 0.5

    search_radius = 10.
    idx1, idx2, d2d, d_ra, d_dec = match_coord([ra], [dec], exp['ra_bore'], exp['dec_bore'], search_radius=search_radius, plot_q=False, keep_all_pairs=True, verbose=False)
    nexp = len(idx2)

    if len(idx2)>120:
        tmp = exp[idx2].copy()
        tmp['score'] = tmp['efftime'] * tmp['weight']
        tmp.sort('score')
        tmp = tmp[:-120]
        expnum_list_to_discard += list(tmp['expnum'])

print('Exposures discarded', len(expnum_list_to_discard), len(expnum_list_to_discard)-len(np.unique(expnum_list_to_discard)))
expnum_list_to_discard = np.unique(expnum_list_to_discard)

# mask = np.in1d(exp['expnum'], expnum_list_to_discard)
# exp = exp[~mask]
mask = np.in1d(ccd['expnum'], expnum_list_to_discard)
ccd = ccd[~mask]

#######################################################################################################

def get_nccd(ccd_index):
    mask = (bricks['RA']>ccd['ra'][ccd_index]-half_width/np.cos(np.radians(bricks['DEC']))) & (bricks['RA']<ccd['ra'][ccd_index]+half_width/np.cos(np.radians(bricks['DEC']))) \
        & (bricks['DEC']>ccd['dec'][ccd_index]-half_height) & (bricks['DEC']<ccd['dec'][ccd_index]+half_height)
    nccd = np.array(mask, dtype=int)
    return nccd


def get_brickscore(ccd_index):
    mask = (bricks['RA']>ccd['ra'][ccd_index]-half_width/np.cos(np.radians(bricks['DEC']))) & (bricks['RA']<ccd['ra'][ccd_index]+half_width/np.cos(np.radians(bricks['DEC']))) \
        & (bricks['DEC']>ccd['dec'][ccd_index]-half_height) & (bricks['DEC']<ccd['dec'][ccd_index]+half_height)
    if np.sum(mask)==0:
        return -np.inf
    if np.all(bricks['nccd'][mask]<nccd_threshold):
        return -np.inf
    brickscore = np.sum(bricks['nccd'][mask]-nccd_threshold)
    return brickscore


ccds_to_keep = np.arange(len(ccd))

with Pool(processes=n_processes) as pool:
    res = np.array(pool.map(get_nccd, ccds_to_keep))
bricks['nccd'] = np.sum(res, axis=0)
# bricks['efftime'] = np.dot(res.T, ccd['efftime'][ccds_to_keep])
# bricks['depth'] = 2.5*np.log10(np.sqrt(bricks['efftime']/des_nominal_efftimes[band]))

counter = 0
while np.sum(bricks['nccd']>nccd_threshold)>nbricks_threshold:
    counter += 1
    print('Iteration {}'.format(counter))

    with Pool(processes=n_processes) as pool:
        brickscore = np.array(pool.map(get_brickscore, ccds_to_keep))

    actual_step_size = np.minimum(np.max(bricks['nccd']-nccd_threshold), step_size)

    ccds_to_delete = np.argsort(brickscore/ccd['efftime'][ccds_to_keep])[-actual_step_size:]
    ccds_to_keep = np.delete(ccds_to_keep, ccds_to_delete)
    print('nccd', len(ccds_to_keep))

    with Pool(processes=n_processes) as pool:
        res = np.array(pool.map(get_nccd, ccds_to_keep))
    bricks['nccd'] = np.sum(res, axis=0)
    # bricks['efftime'] = np.dot(res.T, ccd['efftime'][ccds_to_keep])
    # bricks['depth'] = 2.5*np.log10(np.sqrt(bricks['efftime']/des_nominal_efftimes[band]))
    print('Above-threshold bricks:', np.sum(bricks['nccd']>nccd_threshold))


print(len(ccds_to_keep))

ccd = ccd[ccds_to_keep]
ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/tmp/survey-ccds-{}-{}.fits'.format(field_index, band), overwrite=True)
