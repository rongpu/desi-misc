from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits
from astropy import wcs

from multiprocessing import Pool
import argparse

from scipy.interpolate import RectBivariateSpline
from scipy.ndimage.filters import gaussian_filter
from scipy import stats

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'}
plt.rcParams.update(params)

n_node = 1 # Haswell
n_processess = 10

# parser = argparse.ArgumentParser()
# parser.add_argument('task_id')
# args = parser.parse_args()
# task_id = int(args.task_id)

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

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

fringe_old_dir = '/global/homes/d/djschleg/cosmo/staging/decam/DECam_CP-Fringe'
fringe_new_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/DECam_CP-Fringe'
image_dir = '/global/project/projectdirs/cosmo/staging/'
# surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'
surveyccd_path = '/global/homes/r/rongpu/mydesi/dr9/fringe/misc/survey-ccds-decam-dr9-z-band-only-trim.fits'
blob_dir = '/global/cscratch1/sd/rongpu/fringe/decam_ccd_blob_mask'

# frgscale_dir = '/global/cscratch1/sd/rongpu/fringe/tmp_frgscale/'
frgscale_dir = '/global/cscratch1/sd/rongpu/fringe/frgscale/'

# image_output_dir = '/global/cscratch1/sd/rongpu/fringe/fringe_corrected_image/'
# image_output_dir = '/global/cscratch1/sd/rongpu/fringe/tmp_img/'
image_output_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/fringe_corrected_image'

##############################################################################################################################

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts']
ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
# ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
mask = ccd['ccd_cuts']==0
mask &= ccd['filter']=='z' # include only z-band images
ccd = ccd[mask]
print(len(ccd))
ccd['ccdnum'] = [ccdnamenumdict[ccd['ccdname'][ii].strip()] for ii in range(len(ccd))]

expnum_list = np.unique(ccd['expnum'])
# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)

# # split among the Cori nodes
# expnum_list_split = np.array_split(expnum_list, n_node)
# expnum_list = expnum_list_split[task_id]

# Load old fringe image
fringe_old_dict = {}
for ccdnum in range(1, 63):
    # skip N30 and S7
    if ccdnum==61 or ccdnum==31:
        continue
    fringe_old_path = os.path.join(fringe_old_dir, 'DES17B_20180103_908c062-z-{}_frg.fits'.format(str(ccdnum).zfill(2)))
    fringe_old = fits.getdata(fringe_old_path)
    # remove the edge pixels
    fringe_old = fringe_old[1:4095, 1:2047]
    fringe_old_dict[ccdnum] = fringe_old.copy()

# Load new fringe images
fringe_new_dict = {}
for ccdnum in range(1, 63):
    # skip N30 and S7
    if ccdnum==61 or ccdnum==31:
        continue
    fringe_path = glob.glob(os.path.join(fringe_new_dir, '*CCD{}.fits'.format(str(ccdnum).zfill(2))))[0]
    fringe_img = fits.getdata(fringe_path)
    fringe_img = fringe_img[1:4095, 1:2047]
    fringe_new_dict[ccdnum] = fringe_img.copy()

##############################################################################################################################

def save_image(expnum):

    print('expnum:', expnum)

    # Find an arbitrary CCD in the exposure to get the image filename
    ccd_index = np.where((ccd['expnum']==expnum))[0][0]

    frgscale_path = os.path.join(frgscale_dir, ccd['image_filename'][ccd_index].strip().replace('.fits.fz', '.txt'))
    image_output_path = os.path.join(image_output_dir, ccd['image_filename'][ccd_index].strip())

    if os.path.isfile(image_output_path):
        print(image_output_path+' already exists!')
        return None

    if not os.path.exists(os.path.dirname(image_output_path)):
        try:
            os.makedirs(os.path.dirname(image_output_path))
        except:
            pass

    try:
        fringe_table = Table.read(frgscale_path, format='ascii.commented_header')
        fringe_table.rename_column('slope', 'frgscale')
        mask = fringe_table['ccdnum']!=2
        n_ccd = np.sum(mask)
        frgscale_median = np.median(fringe_table['frgscale'][mask])
    except:
        n_ccd = 0


        #!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%
        print('fringe fitting results not available!')
        return None
        #!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%#!#@$%


    img_fn = os.path.join(image_dir, ccd['image_filename'][ccd_index].strip())
    hdulist = fits.open(img_fn)

    for hdu_index in range(1, len(hdulist)):

        # print(hdu_index)

        ccdname = hdulist[hdu_index].header['EXTNAME'].strip()
        ccdnum = ccdnamenumdict[ccdname]

        img = hdulist[hdu_index].data

        try:
            frgscale_old = (hdulist[hdu_index].header)['FRGSCALE']
        except:
            continue

        # Back out the exisiting fringe correction
        fringe_old = fringe_old_dict[ccdnum]
        img += fringe_old*frgscale_old

        # require at least 10 CCDs for fringe correction
        if n_ccd<10:
            print('HDU{}: Only {} CCDs available; Fringe correction not performed.'.format(hdu_index, n_ccd))
            hdulist[hdu_index].data = img
            continue
        
        mask = fringe_table['ccdnum']==ccdnum
        if np.sum(mask)!=1:
            continue
        table_index = np.where(mask)[0][0]
        frgscale_new = fringe_table['frgscale'][table_index]
        fringe_new = fringe_new_dict[ccdnum]

        # Fring scale shall not be larger than 3 times the median value
        if frgscale_new>3*frgscale_median:
            frgscale_new = 3*frgscale_median

        img -= frgscale_new * fringe_new

        hdulist[hdu_index].data = img

    hdulist.writeto(image_output_path)
    hdulist.close()

    gc.collect()

def main():

    # with Pool(processes=n_processess) as pool:
    #     res = pool.map(save_image, expnum_list)

    for expnum in expnum_list:
        save_image(expnum)

    # save_image(776310)

if __name__=="__main__":
    main()

