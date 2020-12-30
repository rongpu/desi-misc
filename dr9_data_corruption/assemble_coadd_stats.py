from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

from multiprocessing import Pool

n_processes = 32

field = 'south'
output_dir = '/global/cscratch1/sd/rongpu/dr9_tests/bad_coadds'

bricks = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/survey-bricks.fits.gz')
brick_index_all = np.arange(len(bricks))

for band in ['g', 'r', 'z']:
    bricks['npix_obs_'+band] = -99
    bricks['npix_case1_'+band] = -99
    bricks['npix_case2_'+band] = -99
    bricks['npix_case3_'+band] = -99

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
brick_index_all = np.random.choice(brick_index_all, size=len(brick_index_all), replace=False)

brick_index_split = np.array_split(brick_index_all, n_processes)

def assemble_coadd_stats(brick_index_list):

    cat = bricks[brick_index_list]

    for cat_index in range(len(cat)):

        brickname = cat['BRICKNAME'][cat_index]
        output_path = os.path.join(output_dir, '{}'.format(brickname[:3]), 'legacysurvey-{}.txt'.format(brickname))

        if not os.path.isfile(output_path):
            continue

        output = np.loadtxt(output_path).astype(int)
        for index, band in enumerate(['g', 'r', 'z']):
            cat['npix_obs_'+band][cat_index], cat['npix_case1_'+band][cat_index], cat['npix_case2_'+band][cat_index], cat['npix_case3_'+band][cat_index] = \
            output[index]

    mask = cat['npix_obs_r']>0
    cat = cat[mask]

    return cat

def main():

    print('Starting')

    with Pool(processes=n_processes) as pool:
        res = pool.map(assemble_coadd_stats, brick_index_split)

    cat = vstack(res)
    cat.write('/global/cfs/cdirs/desi/users/rongpu/dr9/outliers/coadd_stats.fits')

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

