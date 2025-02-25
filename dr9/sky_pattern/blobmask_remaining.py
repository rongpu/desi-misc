# Create a survey-ccd catalog of the remaining exposures yet to hae blobmask

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'filter', 'mjd_obs', 'ra', 'dec', 'ccdraoff', 'ccddecoff', 'ccd_cuts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))

ccd_all = ccd.copy()

# survey-ccd with unique exposures only
ccd.sort('expnum')
mask = np.concatenate([[True], np.diff(ccd['expnum'])!=0])
ccd = ccd[mask]

exposures_done = np.zeros(len(ccd), dtype=bool)
for ccd_index in range(len(ccd)):
    str_loc = str.find(ccd['image_filename'][ccd_index].strip(), '.fits')
    img_filename_base = ccd['image_filename'][ccd_index].strip()[:str_loc]
    blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
    if os.path.isfile(blob_path):
        exposures_done[ccd_index] = True

print(np.sum(exposures_done), len(ccd), np.sum(exposures_done)/len(ccd))

expnum_remain = ccd['expnum'][~exposures_done]

mask = np.in1d(ccd_all['expnum'], expnum_remain)
ccd_remain = ccd_all[mask]
ccd_remain.write('/global/cscratch1/sd/rongpu/temp/survey-ccds-decam-dr9-remaining.fits')
