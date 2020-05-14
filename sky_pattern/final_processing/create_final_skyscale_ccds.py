from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits

from multiprocessing import Pool

n_processes = 32

skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8.fits')

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
ccd_columns = ['image_hdu', 'expnum', 'ccdname', 'filter', 'ccdskycounts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))

ccd['ccd_id'] = 100*ccd['expnum'] + ccd['image_hdu']

skyscale = Table.read('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_raw.fits')
skyscale.remove_column('ccdname')

ccd_new = join(ccd, skyscale, keys=['image_hdu', 'expnum'], join_type='left')
print(len(ccd))
print(len(ccd_new))

mask = np.copy(ccd_new['ccdskyscale'].mask)
ccd_new = ccd_new.filled(fill_value=0)
ccd_new['ccdskyscale'][mask] = 0
ccd_new['medianskyscale'][mask] = 0
ccd_new['run'][mask] = -1

ccd_new['skyscale'] = 0.

expnum_list = np.unique(ccd_new['expnum'])

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)

def get_ccds(expnum):

    # if expnum%100==0:
    #     print(expnum, len(expnum_list))

    # All CCDs in the exposure
    mask_exp = (ccd_new['expnum']==expnum)
    ccd_exp = (ccd_new[mask_exp]).copy()

    if not expnum in skyrun['expnum']:
        return ccd_exp

    skyrun_index = np.where(skyrun['expnum']==expnum)[0][0]
    run = skyrun['run'][skyrun_index]

    # CCDs with valid skyscales measured
    mask_has_skyscale = (ccd_exp['run']>=0)

    if np.sum(mask_has_skyscale)>=10:
        ccd_exp['skyscale'] = (ccd_exp['medianskyscale'][mask_has_skyscale])[0]
    else:
        ccd_exp['skyscale'] = np.median(ccd_exp['ccdskycounts'])

    # Fill in the run numbers last
    ccd_exp['run'] = run

    return ccd_exp


def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(get_ccds, expnum_list)

    print('pool done!')
    print(len(res))

    ccd_new1 = vstack(res)
    print(len(ccd_new1))

    # Line match to survey-ccds
    ccd_reverse_sort = np.array(ccd['ccd_id']).argsort().argsort()
    ccd_new1.sort('ccd_id')
    ccd_new1 = ccd_new1[ccd_reverse_sort]

    ccd_new.write('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_debug.fits')

    # Remove the skyscale values that won't be used
    ccd_new1.remove_columns(['ccdskyscale', 'medianskyscale', 'ccdskycounts', 'ccd_id'])

    ccd_new1.write('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds.fits')

    print('All done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

