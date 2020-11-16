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


n_processes = 32

# field = 'north'
field = 'south'
coadd_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/{}/coadd'.format(field)

bricks = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/survey-bricks.fits.gz')
brickid_all = np.unique(bricks['BRICKID'])

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
brickid_all = np.random.choice(brickid_all, size=len(brickid_all), replace=False)

# split among the processes
brickid_split = np.array_split(brickid_all, n_processes)


def get_wise_sky_residuals(brickid_list):

    brickid_list = np.sort(brickid_list)

    mask = np.in1d(bricks['BRICKID'], brickid_list)
    residuals = bricks[['BRICKNAME', 'BRICKID', 'RA', 'DEC']][mask].copy()
    residuals.sort('BRICKID')

    if not np.all(residuals['BRICKID']==brickid_list):
        raise ValueError('BRICKID no do match!')

    columns = ['w1_median', 'w1_mean', 'w1_nmad', 'w1_10th_percentile', 'w1_90th_percentile', 'w2_median', 'w2_mean', 'w2_nmad', 'w2_10th_percentile', 'w2_90th_percentile']
    for col in columns:
        residuals[col] = -99.

    for brickid in brickid_list:

        index = np.where(residuals['BRICKID']==brickid)[0][0]
        brickname = residuals['BRICKNAME'][index]
        data_dir = os.path.join(coadd_dir, '{}/{}'.format(brickname[:3], brickname))

        model_path = os.path.join(data_dir, 'legacysurvey-{}-model-{}.fits.fz'.format(brickname, 'W1'))
        if not os.path.isfile(model_path):
            continue

        for band in ['W1', 'W2']:

            # img_path = glob.glob(os.path.join(data_dir, '*-image-{}.fits.fz'.format(band)))[0]
            # model_path = glob.glob(os.path.join(data_dir, '*-model-{}.fits.fz'.format(band)))[0]
            img_path = os.path.join(data_dir, 'legacysurvey-{}-image-{}.fits.fz'.format(brickname, band))
            model_path = os.path.join(data_dir, 'legacysurvey-{}-model-{}.fits.fz'.format(brickname, band))

            img = fitsio.read(img_path)
            model = fitsio.read(model_path)
            # print(img.shape)

            npix_trim = 16
            mask = np.zeros_like(img, dtype=bool)
            mask[npix_trim:-npix_trim, npix_trim:-npix_trim] = True
            img[~mask] = np.nan

            img_flat = (img-model).flatten()
            img_flat = img_flat[np.isfinite(img_flat)]

            residuals[band.lower()+'_median'][index] = np.median(img_flat)
            residuals[band.lower()+'_mean'][index] = np.mean(img_flat)
            residuals[band.lower()+'_nmad'][index] = nmad(img_flat)
            residuals[band.lower()+'_10th_percentile'][index] = np.percentile(img_flat, 10)
            residuals[band.lower()+'_90th_percentile'][index] = np.percentile(img_flat, 90)

    return residuals


def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(get_wise_sky_residuals, brickid_split)

    print('All done!!!!!!!!!!!!!!!')

    residuals = vstack(res)
    residuals.sort('BRICKID')

    residuals.write('/global/u2/r/rongpu/data/survey_bricks_wise_sky_residuals-dr9_{}.fits'.format(field))


if __name__=="__main__":
    main()
