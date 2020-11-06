# This code is too slow to finish running in 4 hours on an interactive node

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

from multiprocessing import Pool


def nmad(x): return 1.4826 * np.median(np.abs(x-np.median(x)))


n_processess = 32

field = 'north'
# field = 'south'
coadd_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/{}/coadd'.format(field)
metrics_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/{}/metrics'.format(field)

bits = [0, 1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 13]

bricks = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/survey-bricks.fits.gz')
brickid_all = np.unique(bricks['BRICKID'])

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
brickid_all = np.random.choice(brickid_all, size=len(brickid_all), replace=False)

# split among the processes
brickid_split = np.array_split(brickid_all, n_processess)


def get_optical_sky_residuals(brickid_list):

    brickid_list = np.sort(brickid_list)

    mask = np.in1d(bricks['BRICKID'], brickid_list)
    residuals = bricks[['BRICKNAME', 'BRICKID', 'RA', 'DEC']][mask].copy()
    residuals.sort('BRICKID')

    if not np.all(residuals['BRICKID']==brickid_list):
        raise ValueError('BRICKID no do match!')

    columns = ['frac', 'blobmask_frac', 'g_median', 'g_mean', 'g_nmad', 'g_blobmask_median', 'g_blobmask_mean', 'g_blobmask_nmad', 'r_median', 'r_mean', 'r_nmad', 'r_blobmask_median', 'r_blobmask_mean', 'r_blobmask_nmad', 'z_median', 'z_mean', 'z_nmad', 'z_blobmask_median', 'z_blobmask_mean', 'z_blobmask_nmad']
    for col in columns:
        residuals[col] = -99.

    for brickid in brickid_list:

        index = np.where(residuals['BRICKID']==brickid)[0][0]
        brickname = residuals['BRICKNAME'][index]
        data_dir = os.path.join(coadd_dir, '{}/{}'.format(brickname[:3], brickname))

        blobmodel_path = os.path.join(data_dir, 'legacysurvey-{}-blobmodel-{}.fits.fz'.format(brickname, 'g'))
        if not os.path.isfile(blobmodel_path):
            continue

        maskbits_path = os.path.join(data_dir, 'legacysurvey-{}-maskbits.fits.fz'.format(brickname))
        blob_path = glob.glob(os.path.join(metrics_dir, '{}/blobs-{}.fits.gz'.format(brickname[:3], brickname)))[0]
        if (not os.path.isfile(maskbits_path)) or (not os.path.isfile(blob_path)):
            continue

        three_band_coverage = True
        for band in ['g', 'r', 'z']:
            blobmodel_path = os.path.join(data_dir, 'legacysurvey-{}-blobmodel-{}.fits.fz'.format(brickname, band))
            if not os.path.isfile(blobmodel_path):
                three_band_coverage = False
        if not three_band_coverage:
            continue

        maskbits = fitsio.read(maskbits_path)
        blob = fitsio.read(blob_path)
        blob = blob!=-1  # True means within a blob
        # print(img.shape)

        mask_clean = np.ones_like(maskbits, dtype=bool)
        for bit in bits:
            mask_clean &= (maskbits & 2**bit)==0
        # print(np.sum(~mask_clean)/np.product(mask_clean.shape))

        residuals['frac'][index] = np.sum(mask_clean)/np.product(mask_clean.shape)
        residuals['blobmask_frac'][index] = np.sum(mask_clean & (~blob))/np.product(mask_clean.shape)

        # print(residuals['frac'][index])
        # print(residuals['blobmask_frac'][index])

        if residuals['frac'][index]==0:
            continue

        for band in ['g', 'r', 'z']:

            # img_path = glob.glob(os.path.join(data_dir, '*-image-{}.fits.fz'.format(band)))[0]
            # blobmodel_path = glob.glob(os.path.join(data_dir, '*-blobmodel-{}.fits.fz'.format(band)))[0]
            img_path = os.path.join(data_dir, 'legacysurvey-{}-image-{}.fits.fz'.format(brickname, band))
            blobmodel_path = os.path.join(data_dir, 'legacysurvey-{}-blobmodel-{}.fits.fz'.format(brickname, band))

            try:
                img = fitsio.read(img_path)
            except OSError:
                print('File not found: {}'.format(img_path))
                continue
            try:
                blobmodel = fitsio.read(blobmodel_path)
            except OSError:
                print('File not found: {}'.format(blobmodel_path))
                continue

            img[~mask_clean] = np.nan
            img_flat = (img-blobmodel).flatten()
            img_flat = img_flat[np.isfinite(img_flat)]

            residuals[band.lower()+'_median'][index] = np.median(img_flat)
            residuals[band.lower()+'_mean'][index] = np.mean(img_flat)
            residuals[band.lower()+'_nmad'][index] = nmad(img_flat)

            if residuals['blobmask_frac'][index]==0:
                continue

            img[blob] = np.nan
            img_flat = (img-blobmodel).flatten()
            img_flat = img_flat[np.isfinite(img_flat)]

            residuals[band.lower()+'_blobmask_median'][index] = np.median(img_flat)
            residuals[band.lower()+'_blobmask_mean'][index] = np.mean(img_flat)
            residuals[band.lower()+'_blobmask_nmad'][index] = nmad(img_flat)

    return residuals


def main():

    with Pool(processes=n_processess) as pool:
        res = pool.map(get_optical_sky_residuals, brickid_split)

    print('All done!!!!!!!!!!!!!!!')

    residuals = vstack(res)
    residuals.sort('BRICKID')

    residuals.write('/global/u2/r/rongpu/data/survey_bricks_optical_sky_residuals-dr9_{}.fits'.format(field))


if __name__=="__main__":
    main()
