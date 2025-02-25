from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from multiprocessing import Pool
import argparse
from pathlib import Path

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'}
plt.rcParams.update(params)

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

fringe_norm_dict = {1: 1.074, 2: 1.500, 3: 1.010, 4: 1.120, 5: 1.080, 6: 1.130, 7: 1.023, 
8: 1.016, 9: 1.018, 10: 0.943, 11: 1.000, 12: 1.000, 13: 1.002, 14: 1.030, 15: 1.040, 
16: 1.244, 17: 1.019, 18: 1.082, 19: 0.951, 20: 0.979, 21: 0.958, 22: 0.948, 23: 0.971, 
24: 1.113, 25: 0.950, 26: 0.956, 27: 0.963, 28: 0.998, 29: 0.960, 30: 0.986, 32: 1.058, 
33: 0.943, 34: 0.992, 35: 0.974, 36: 0.949, 37: 1.015, 38: 1.031, 39: 1.042, 40: 1.028, 
41: 0.944, 42: 0.938, 43: 0.978, 44: 1.112, 45: 1.039, 46: 0.916, 47: 0.997, 48: 0.883, 
49: 1.019, 50: 1.017, 51: 1.016, 52: 1.009, 53: 1.007, 54: 0.950, 55: 0.913, 56: 1.084, 
57: 0.785, 58: 0.946, 59: 0.981, 60: 1.049, 62: 1.024}

fringe_old_dir = '/global/homes/d/djschleg/cosmo/staging/decam/DECam_CP-Fringe'
# Here using the normalized fringe template, unlike in fringe_template_fitting.py:
fringe_new_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/DECam_CP-Fringe-Normed'
image_dir = '/global/project/projectdirs/cosmo/staging/'
surveyccd_path = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/temp/survey-ccds-decam-dr9-z-band-only-trim-all-dr9.fits'
# surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'
# surveyccd_path = '/global/homes/r/rongpu/mydesi/dr9/fringe/misc/survey-ccds-decam-dr9-z-band-only-trim.fits'
blob_dir = '/global/cscratch1/sd/rongpu/fringe/decam_ccd_blob_mask'

frgscale_dir = '/global/cscratch1/sd/rongpu/fringe/frgscale/'

frgscale_output_dir = '/global/cscratch1/sd/rongpu/fringe/frgscale_applied_20200302/'
image_output_dir = '/global/cscratch1/sd/rongpu/fringe/fringe_corrected_image_20200302/'

status_dir = '/global/homes/r/rongpu/temp/status'

# image_output_dir = '/global/cscratch1/sd/rongpu/fringe/tmp_img/'
# image_output_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/fringe_corrected_image'

##############################################################################################################################

# Load CCD list
# ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts', 'mjd_obs']
# ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
# mask = ccd['ccd_cuts']==0
# mask &= ccd['filter']=='z' # include only z-band images
mask = ccd['filter']=='z' # include only z-band images
ccd = ccd[mask]
print(len(ccd))
ccd['ccdnum'] = [ccdnamenumdict[ccd['ccdname'][ii].strip()] for ii in range(len(ccd))]
# expnum_list = np.unique(ccd['expnum'])

# expnum_list_delete = [230089, 232685, 243126, 251361, 253121, 253471, 263067, 277619, 280310, 281616, 283436, 298157, 302268, 302671, 302757, 314055, 314063, 314157, 314318, 318938, 318951, 323273, 335357, 339925, 339926, 339927, 339931, 339943, 339946, 348008, 355320, 358514, 359348, 359693, 367541, 371338, 374599, 374874, 376734, 378903, 382167, 382230, 384064, 385662, 385664, 385665, 385701, 393598, 396652, 404845, 405884, 405988, 417101, 417156, 421993, 422531, 422541, 422554, 423357, 423714, 426050, 429358, 429360, 429836, 429842, 429850, 430285, 430353, 448695, 448699, 448700, 458365, 477253, 477306, 492438, 494960, 494982, 497199, 497206, 502448, 504579, 511344, 513311, 515687, 515688, 515695, 515697, 515701, 519029, 520558, 521770, 521887, 522223, 524799, 529077, 529832, 529833, 530184, 533702, 533708, 534675, 534683, 564019, 572501, 591508, 593514, 609588, 616006, 618264, 618286, 618294, 618296, 618297, 618299, 618308, 618312, 618313, 618647, 618650, 618654, 618656, 618660, 618662, 618664, 618678, 618682, 619026, 619027, 619028, 619040, 619673, 619675, 619676, 619680, 619698, 619703, 619720, 672431, 680057, 684280, 693008, 695115, 703207, 719612, 723655, 774115, 775406, 775936, 783884, 784332, 808667, 865273, 872515, 872523, 872530, 872536, 872537, 872563, 872571, 872575, 872577, 873030, 873034, 873053, 873058, 873061, 873066, 873074, 873086, 873092, 873095, 873101, 873104, 873106, 873114, 873117, 873124, 873126, 873131, 873144, 873145, 873152, 873155, 873159, 873183, 873442, 873446, 873448, 873457, 873458, 873465, 873466, 873471, 873472, 873474, 873478, 873488, 873494, 873502, 873506, 873507, 873516, 873519, 873523, 873524, 873525, 873538, 873539, 873543]

expnum_list_delete = [254222, 256596, 267230, 346644, 347144, 347976, 359694, 359753, 376761, 382216, 386384, 386745, 405932, 425653, 430734, 431718, 472398, 486322, 488811, 504603, 506061, 506072, 513357, 513382, 515962, 521345, 522187, 548231, 553437, 563055, 564477, 573146, 573375, 581506, 588877, 589754, 591898, 592293, 593155, 593476, 601789, 603129, 618664, 625231, 625664, 631028, 642215, 647421, 648460, 672368, 681941, 681956, 690655, 693115, 695099, 704147, 705183, 718885, 719415, 725182, 746182, 763979, 768048, 768493, 769228, 874277]

print(len(expnum_list_delete))

counter = 0

for expnum in expnum_list_delete:

    ccd_index = np.where((ccd['expnum']==expnum))[0][0]
    img_fn_write = os.path.join(image_output_dir, ccd['image_filename'][ccd_index].strip())

    if os.path.isfile(img_fn_write):
        print(img_fn_write)
        os.remove(img_fn_write)
        counter += 1
print(counter)
