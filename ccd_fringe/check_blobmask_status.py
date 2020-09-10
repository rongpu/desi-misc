from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
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

surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-trim.fits'
# surveyccd_path_dr8 = '/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr8-trim.fits'
blob_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'

# Load CCD list
ccd_columns = ['image_filename', 'expnum']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
print(len(ccd))

# ccd = fitsio.read(surveyccd_path)
# mask = ccd['ccd_cuts']==0
# band = 'z'
# mask = ccd['filter']==band
# ccd = ccd[mask]
# print(len(ccd))

status_arr = np.zeros(len(ccd), dtype=bool)

for index, image_filename in enumerate(ccd['image_filename']):

    image_filename = image_filename.strip()
    str_loc = str.find(image_filename, '.fits')
    img_filename_base = image_filename[:str_loc]
    blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
    if os.path.isfile(blob_path):
        status_arr[index] = True

print(np.sum(status_arr), np.sum(~status_arr), np.sum(status_arr)/len(status_arr))

