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

ccdnum_list = [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
       18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
       35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
       52, 53, 54, 55, 56, 57, 58, 59, 60, 62]

ccdnamenumdict = {'S1': 25, 'S2': 26, 'S3': 27, 'S4':28,
                  'S5': 29, 'S6': 30, 'S7': 31,
                  'S8': 19, 'S9': 20, 'S10': 21, 'S11': 22, 'S12': 23,
                  'S13': 24,
                  'S14': 13, 'S15': 14, 'S16': 15, 'S17': 16, 'S18': 17,
                  'S19': 18,
                  'S20': 8, 'S21': 9, 'S22': 10, 'S23': 11, 'S24': 12,
                  'S25': 4, 'S26': 5, 'S27': 6, 'S28': 7,
                  'S29': 1, 'S30': 2, 'S31': 3,
                  'N1': 32, 'N2': 33, 'N3': 34, 'N4': 35,
                  'N5': 36, 'N6': 37, 'N7': 38,
                  'N8': 39, 'N9': 40, 'N10': 41, 'N11': 42, 'N12': 43,
                  'N13': 44,
                  'N14': 45, 'N15': 46, 'N16': 47, 'N17': 48, 'N18': 49,
                  'N19': 50,
                  'N20': 51, 'N21': 52, 'N22': 53, 'N23': 54, 'N24': 55,
                  'N25': 56, 'N26': 57, 'N27': 58, 'N28': 59,
                  'N29': 60, 'N30': 61, 'N31': 62,
                  }
ccdnamenumdict_inv = {aa: bb for bb, aa in ccdnamenumdict.items()}


n_processess = 32

skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8.fits')
run_list = np.unique(skyrun['run'])

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
run_list = np.random.choice(run_list, size=len(run_list), replace=False)

def patch_s7(run):
    
    mask = skyrun['run']==run
    band = skyrun['filter'][mask][0]

    print('band: {}, run: {}'.format(band, run))

    img_all_path = os.path.join('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates_v2/sky_template_{}_{}.fits.fz'.format(band, run))
    img_s7_path = os.path.join('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates_s7/sky_template_{}_{}.fits.fz'.format(band, run))
    img_all_new_path = os.path.join('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates_v2_patched/sky_template_{}_{}.fits.fz'.format(band, run))

    if os.path.isfile(img_all_new_path):
        print(img_all_new_path, 'already exist!')
        return None

    # # The file should be at least 0.3 hours old to ensure it's not being written
    # time_modified = os.path.getmtime(img_all_path)
    # if (time.time() - time_modified)/3600 < 0.3:
    #     print('file too new; skipping')
    #     return None

    hdul_w = fitsio.FITS(img_all_new_path, mode='rw', clobber=True)
    hdul_w.write(data=None) # first HDU is empty
    
    for ccdnum in ccdnum_list:

        ccdname = ccdnamenumdict_inv[ccdnum]

        if ccdname=='S7':
            img = fitsio.read(img_s7_path, ext=ccdname)
        else:
            try:
                img = fitsio.read(img_all_path, ext=ccdname)
            except:
                print('ERROR!!!', 'band: {}, run: {}'.format(band, run), '- Could not read {}'.format(ccdname))
                continue

        hdul_w.write(data=img, extname=ccdname, compress='rice')

    hdul_w.close()


def main():

    with Pool(processes=n_processess) as pool:
        res = pool.map(patch_s7, run_list)

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

