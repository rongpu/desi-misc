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

image_dir = '/global/project/projectdirs/cosmo/staging'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'

skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8.fits')

ccd_columns = ['image_filename', 'expnum', 'ccdname', 'filter', 'ccdskycounts', 'plver']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
ccd['ccdname'] = np.char.strip(ccd['ccdname']) # strip spaces
ccd['plver'] = np.char.strip(ccd['plver']) # strip spaces

mask = ccd['plver']=='V4.9'
print(np.sum(mask))
ccd = ccd[mask]

# ccd['ccd_id'] = 100*ccd['expnum'] + ccd['image_hdu']

skyscale = Table.read('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_raw_cp_v4.9.fits')
skyscale.remove_column('image_hdu')

ccd_new = join(ccd, skyscale, keys=['ccdname', 'expnum'], join_type='left')
print(len(ccd))
print(len(ccd_new))

mask = np.copy(ccd_new['ccdskyscale'].mask)
ccd_new = ccd_new.filled(fill_value=0)
ccd_new['ccdskyscale'][mask] = 0
ccd_new['medianskyscale'][mask] = 0
ccd_new['run'][mask] = -1

ccd_new['skyscale'] = 0.

expnum_list = np.unique(ccd_new['expnum'])

# # shuffle
# np.random.seed(123)
# # DO NOT USE NP.RANDOM.SHUFFLE
# expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)

def get_ccds(expnum):

    # if expnum%100==0:
    #     print(expnum, len(expnum_list))

    # All CCDs in the exposure
    mask_exp = (ccd_new['expnum']==expnum)
    ccd_exp = (ccd_new[mask_exp]).copy()

    if not expnum in skyrun['expnum']:
        ccd_exp['PLPROCID'] = ''
        return ccd_exp

    skyrun_index = np.where(skyrun['expnum']==expnum)[0][0]
    run = skyrun['run'][skyrun_index]
    image_filename = ccd_exp['image_filename'][0].strip()

    # CCDs with valid skyscales measured
    mask_has_skyscale = (ccd_exp['run']>=0)

    if np.sum(mask_has_skyscale)>=10:
        ccd_exp['skyscale'] = (ccd_exp['medianskyscale'][mask_has_skyscale])[0]
    else:
        ccd_exp['skyscale'] = np.median(ccd_exp['ccdskycounts'])

    # Fill in the run numbers last
    ccd_exp['run'] = run

    # Add PLPROCID
    img_path = os.path.join(image_dir, image_filename)
    with fits.open(img_path) as f:
        ccd_exp['PLPROCID'] = f[0].header['PLPROCID']

    return ccd_exp


def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(get_ccds, expnum_list)

    print('pool done!')
    print(len(res))

    ccd_new1 = vstack(res)
    print(len(ccd_new1))

    # # Line match to survey-ccds
    # ccd_reverse_sort = np.array(ccd['ccd_id']).argsort().argsort()
    # ccd_new1.sort('ccd_id')
    # ccd_new1 = ccd_new1[ccd_reverse_sort]

    ccd_new1.write('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_debug_cp_v4.9.fits')

    # Remove the skyscale values that won't be used
    ccd_new1.remove_columns(['image_filename', 'plver', 'ccdskyscale', 'medianskyscale', 'ccdskycounts'])

    ccd_new1['SKYTMPLV'] = 2

    ccd_new1.write('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_cp_v4.9.fits')

    print('All done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

