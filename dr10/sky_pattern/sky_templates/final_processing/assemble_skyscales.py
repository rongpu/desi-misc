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

skyscale_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_scales_tmp/'

skyrun = Table(fitsio.read('/global/cfs/cdirs/desi/users/schlafly/decals/skyrunsdr10-v3.fits'))

mask = skyrun['filter']=='z'
skyrun = skyrun[mask]
print('skyrun', len(skyrun))

expnum_list = skyrun['expnum'].copy()


def assemble_skyscales(expnum):

    # Get info
    skyrun_index = np.where(skyrun['expnum']==expnum)[0][0]
    band = skyrun['filter'][skyrun_index]
    run = skyrun['run'][skyrun_index]
    image_filename = skyrun['image_filename'][skyrun_index].strip()

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
    skyscale.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_scales/skyscales_ccds_raw.fits')

    print('All done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()
