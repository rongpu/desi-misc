from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from multiprocessing import Pool

n_processes = 32

image_dir = '/global/cfs/cdirs/cosmo/staging/'

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
# surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9-garage/reorg/decam/survey-ccds-decam-dr8-newlocs2.fits.gz'

ccd = Table(fitsio.read(surveyccd_path, columns=['expnum', 'image_filename', 'filter']))
print(len(ccd))

# Only keep unique exposures
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]
print(len(ccd))

def get_pupil_params(index):

    image_path = os.path.join(image_dir, ccd['image_filename'][index].strip())
    with fitsio.FITS(image_path) as f:
        try:
            header = f[0].read_header()
            if 'PUPILSKY' in header:
                pupilsky = header['PUPILSKY']
            else:
                pupilsky = np.nan
            if 'PUPILMAX' in header:
                pupilmax = header['PUPILMAX']
            else:
                pupilmax = np.nan
            if 'PUPILAMP' in header:
                pupilamp = header['PUPILAMP']
            else:
                pupilamp = np.nan
        except:
            pupilsky = np.nan
            pupilmax = np.nan
            pupilamp = np.nan

    return pupilsky, pupilmax, pupilamp

def main():

    # ccd['PUPILSKY'] = 0.
    # ccd['PUPILMAX'] = 0.
    # ccd['PUPILAMP'] = 0.

    # for index in np.arange(len(ccd)):

    #     if index%100==0:
    #         print(index, '/', len(ccd))

    #     res = get_pupil_params(index)
    #     ccd['PUPILSKY'][index] = res[0]
    #     ccd['PUPILMAX'][index] = res[1]
    #     ccd['PUPILAMP'][index] = res[2]

    # ccd.write('/global/cscratch1/sd/rongpu/dr9dev/pupil_pattern/survey-ccds-decam-dr8-pupil-params.fits')

    # print('Done!!!!!!!!!!!!!!!!!!!!!')

    with Pool(processes=n_processes) as pool:
        res = pool.map(get_pupil_params, np.arange(len(ccd)))

    res = np.array(res)
    ccd['PUPILSKY'] = res[:, 0]
    ccd['PUPILMAX'] = res[:, 1]
    ccd['PUPILAMP'] = res[:, 2]

    ccd.write('/global/cscratch1/sd/rongpu/dr9dev/pupil_pattern/survey-ccds-decam-dr9-pupil-params.fits')

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

