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

skyscale_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales_cp_v4.9/'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'

# ccd_columns = ['image_hdu', 'expnum', 'ccdskycounts']
# ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))

skyrun = Table.read('/global/cfs/cdirs/desi/users/rongpu/dr9/sky_pattern/skyrunsgoodcountexpnumv48dr8.fits')

ccd_columns = ['image_filename', 'expnum', 'plver']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
ccd['plver'] = np.char.strip(ccd['plver']) # strip spaces
mask = ccd['plver']=='V4.9'
ccd = ccd[mask]
print(len(ccd))

expnum_list = np.unique(ccd['expnum'])

def assemble_skyscales(expnum):
    
    # Get info
    skyrun_index = np.where(skyrun['expnum']==expnum)[0][0]
    band = skyrun['filter'][skyrun_index]
    run = skyrun['run'][skyrun_index]

    ccd_index = np.where((ccd['expnum']==expnum))[0][0]
    image_filename = ccd['image_filename'][ccd_index].strip()

    skyscale_path = os.path.join(skyscale_dir, image_filename.replace('.fits.fz', '-skyscale.txt'))
    skyscale = Table.read(skyscale_path, format='ascii.commented_header')

    if len(skyscale)==0:
        return None
    else:
        skyscale['expnum'] = expnum
        skyscale['run'] = run
        return skyscale

def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(assemble_skyscales, expnum_list)

    print('pool done!')
    print(len(res))

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    print('pop done!')
    print(len(res))

    skyscale = vstack(res)
    skyscale.write('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/skyscales_ccds_raw_cp_v4.9.fits')

    print('All done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()
