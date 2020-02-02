# Find exposures left to process

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
sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import search_around

from multiprocessing import Pool
import argparse

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'}
plt.rcParams.update(params)

image_dir = '/global/project/projectdirs/cosmo/staging/'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'
# output_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
# output_dir = '/global/homes/r/rongpu/data/decam_ccd_blob_mask'
output_dir = '/global/cscratch1/sd/rongpu/fringe/decam_ccd_blob_mask'

# Load brick list
bricks = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/survey-bricks.fits.gz')

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'filter', 'mjd_obs', 'ra', 'dec', 'ccdraoff', 'ccddecoff', 'ccd_cuts']
ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
# ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
mask = ccd['ccd_cuts']==0
mask &= ccd['filter']=='z' # include only z-band images
ccd = ccd[mask]
print(len(ccd), 'CCDs')

print('Total nubmer of exposures:', len(np.unique(ccd['expnum'])))

# Find the CCDs whose blobmask files do not yet exist
ccd_mask = np.zeros(len(ccd), dtype=bool) # True if exist
for ccd_index in range(len(ccd)):
    str_loc = str.find(ccd['image_filename'][ccd_index].strip(), '.fits')
    img_filename_base = ccd['image_filename'][ccd_index].strip()[:str_loc]
    blob_path = os.path.join(output_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
    if os.path.isfile(blob_path):
        ccd_mask[ccd_index] = True
print(np.sum(ccd_mask)/len(ccd_mask))
ccd = ccd[~ccd_mask]
print(len(ccd))
expnum_list = np.unique(ccd['expnum'])
print('Nubmer of exposures left to process:', len(expnum_list))
print()

ccd.write('/global/u2/r/rongpu/temp/blobmask_ccd.fits')
