# Normalize the fringe templates so that each exposure the same fringe scale for all its CCDs

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

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
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'
# surveyccd_path = '/global/u2/r/rongpu/temp/survey-ccds-decam-dr9-z-band-only-trim.fits'
blob_dir = '/global/cscratch1/sd/rongpu/fringe/decam_ccd_blob_mask'

frgscale_output_dir = '/global/cscratch1/sd/rongpu/fringe/frgscale/'

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts']
ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
# ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
mask = ccd['ccd_cuts']==0
mask &= ccd['filter']=='z' # include only z-band images
ccd = ccd[mask]
# print(len(ccd))
ccd['ccdnum'] = [ccdnamenumdict[ccd['ccdname'][ii].strip()] for ii in range(len(ccd))]

expnum_list = np.unique(ccd['expnum'])

# Find the exposures that have the fringe scales computed
expnum_list_done = np.zeros(len(expnum_list), dtype=bool)
for index, expnum in enumerate(expnum_list):
    # Find an arbitrary CCD in the exposure to get the image filename
    ccd_index = np.where((ccd['expnum']==expnum))[0][0]
    frgscale_output_path = os.path.join(frgscale_output_dir, ccd['image_filename'][ccd_index].strip().replace('.fits.fz', '.txt'))
    # Check that the file exist and not empty
    if os.path.isfile(frgscale_output_path) and (os.stat(frgscale_output_path).st_size != 0):
        expnum_list_done[index] = True
print(np.sum(expnum_list_done), np.sum(~expnum_list_done), np.sum(expnum_list_done)/len(expnum_list_done))

expnum_list = expnum_list[expnum_list_done]
print(len(expnum_list))

# Randomly select 1000 exposures
np.random.seed(321)
expnum_list = np.random.choice(expnum_list, size=1000, replace=False)
expnum_list.sort()

fringe_all = []

for expnum in expnum_list:

    # Find an arbitrary CCD in the exposure to get the image filename
    ccd_index = np.where((ccd['expnum']==expnum))[0][0]
    img_fn = os.path.join(image_dir, ccd['image_filename'][ccd_index].strip())
    
    frgscale_output_path = os.path.join(frgscale_output_dir, ccd['image_filename'][ccd_index].strip().replace('.fits.fz', '.txt'))
    
    fringe = Table.read(frgscale_output_path, format='ascii.commented_header')
    fringe.rename_column('slope', 'frgscale')
    
    # Remove outlier CCDs from computing median fringe scale
    ccdnum_exclude = [2, 16, 57] # S30, S17, N26
    mask = ~np.in1d(fringe['ccdnum'], ccdnum_exclude)
    
    frgscale_median = np.median(fringe['frgscale'][mask])
    fringe['relative_frg'] = fringe['frgscale'] / frgscale_median
    
    fringe_all.append(fringe)
    
fringe_stack = vstack(fringe_all)
print(len(fringe_stack))

# The normalized templates are obtained by multiplying the current templates with these factors
fringe_template_multiply = {}
print('ccdnum  ccdname norm_factor std')
for ccdnum in np.unique(fringe_stack['ccdnum']):
    mask = fringe_stack['ccdnum']==ccdnum
    norm_factor = np.median(fringe_stack['relative_frg'][mask])
    fringe_template_multiply[ccdnum] = norm_factor
    norm_factor_std = np.std(fringe_stack['relative_frg'][mask])
    print('{:3}     {:3}     {:5.3f}    {:5.3f}'.format(ccdnum, ccdnamenumdict_inv[ccdnum], norm_factor, norm_factor_std))

# print('# ccdnum  norm_factor')
# for ii, jj in fringe_template_multiply.items():
#     print('{:8} {:7.3f}'.format(ii, jj))

# Apply the normalization

fringe_normed_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/DECam_CP-Fringe-Normed'

for ccdnum in np.unique(fringe_stack['ccdnum']):
    ccdname = ccdnamenumdict_inv[ccdnum]
    fringe_path = os.path.join(fringe_new_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
    with fits.open(fringe_path) as hdul:
        hdul[0].data = hdul[0].data * fringe_template_multiply[ccdnum]
        fringe_normed_path = os.path.join(fringe_normed_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
        hdul.writeto(fringe_normed_path)

# # Check that the new templates are indeed different from the old ones
# for ccdnum in [2]:
#     ccdname = ccdnamenumdict_inv[ccdnum]
#     fringe_path = os.path.join(fringe_new_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
#     fringe_normed_path = os.path.join(fringe_normed_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
#     with fits.open(fringe_path) as hdul:
#         data = hdul[0].data
#         plt.imshow(data, vmin=data.min(), vmax=data.max())
#         plt.show()
#     with fits.open(fringe_normed_path) as hdul:
#         data1 = hdul[0].data
#         plt.imshow(data1, vmin=data.min(), vmax=data.max())
#         plt.show()
#         