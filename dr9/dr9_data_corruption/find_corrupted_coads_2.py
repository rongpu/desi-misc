from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

from multiprocessing import Pool

n_processes = 32
n_node = 2
node_index = 1

field = 'south'
coadd_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/{}/coadd'.format(field)

output_dir = '/global/cscratch1/sd/rongpu/dr9_tests/bad_coadds'

bits_to_mask = [0, 2, 3, 4, 5, 6, 7]

bricks = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/survey-bricks.fits.gz')
brick_index_all = np.arange(len(bricks))

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
brick_index_all = np.random.choice(brick_index_all, size=len(brick_index_all), replace=False)

# split among the nodes
brickid_split = np.array_split(brick_index_all, n_node)
brick_index_all = brickid_split[node_index]

def find_bad_coadd(brick_index):

    brickname = bricks['BRICKNAME'][brick_index]
    data_dir = os.path.join(coadd_dir, '{}/{}'.format(brickname[:3], brickname))
    output_path = os.path.join(output_dir, '{}'.format(brickname[:3]), 'legacysurvey-{}.txt'.format(brickname))

    if os.path.isfile(output_path):
        return None

    for band in ['g', 'r', 'z']:
        img_path = os.path.join(data_dir, 'legacysurvey-{}-image-{}.fits.fz'.format(brickname, band))
        if not os.path.isfile(img_path):
            return None

    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))

    output = []

    for band in ['g', 'r', 'z']:

        img = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/south/coadd/{}/{}/legacysurvey-{}-image-{}.fits.fz'.format(brickname[:3], brickname, brickname, band))
        ivar = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/south/coadd/{}/{}/legacysurvey-{}-invvar-{}.fits.fz'.format(brickname[:3], brickname, brickname, band))
        nexp = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/south/coadd/{}/{}/legacysurvey-{}-nexp-{}.fits.fz'.format(brickname[:3], brickname, brickname, band))
        maskbits = fitsio.read('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/south/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(brickname[:3], brickname, brickname))

        npix_obs = np.sum(nexp>0)

        mask_clean = np.ones(img.shape, dtype=bool)
        for bit in bits_to_mask:
            mask_clean &= (maskbits & 2**bit)==0

        # Case 1
        mask = (img==0) & (nexp>0) & (ivar>0)
        mask &= mask_clean
        npix_case1 = np.sum(mask)

        # Case 2
        mask = (np.abs(img)<1e-4) & (nexp>0) & (ivar>0)
        mask &= mask_clean
        # downsize to remove isolated pixels
        binsize = 3
        trim_size_x = mask.shape[1] % binsize
        trim_size_y = mask.shape[0] % binsize
        mask = mask[:(mask.shape[0]-trim_size_y), :(mask.shape[1]-trim_size_x)]
        # to ignore NAN values, use np.nanmean
        mask = np.median(np.median(mask.reshape((mask.shape[0]//binsize, binsize, mask.shape[1]//binsize,-1)), axis=3), axis=1)
        mask = mask>0
        npix_case2 = np.sum(mask)

        # Case 3
        mask = (nexp==0) & (ivar==0) & (img!=0)
        mask &= mask_clean
        npix_case3 = np.sum(mask)

        output.append([npix_obs, npix_case1, npix_case2, npix_case3])

    np.savetxt(output_path, output, fmt='%i')
    # print(output_path)

    return None

def main():

    print('Starting')

    with Pool(processes=n_processes) as pool:
        pool.map(find_bad_coadd, brick_index_all)

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

