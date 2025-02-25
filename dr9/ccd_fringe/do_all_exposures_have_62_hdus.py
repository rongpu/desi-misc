# The answer is no

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

fringe_dir = '/global/homes/d/djschleg/cosmo/staging/decam/DECam_CP-Fringe'
image_dir = '/global/project/projectdirs/cosmo/staging/'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'mjd_obs', 'ra', 'dec', 'skyrms', 'ccdraoff', 'ccddecoff', 'ccd_cuts']
ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
# ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
mask = ccd['ccd_cuts']==0
mask &= ccd['filter']=='z' # include only z-band images
ccd = ccd[mask]
print(len(ccd))

image_filename_unique = np.unique(ccd['image_filename'])

ccd_count = 0
for image_filename in image_filename_unique:
    img_fn = os.path.join(image_dir, image_filename.strip())
    img_hdulist = fits.open(img_fn)
    if len(img_hdulist)!=62:
        print('{} has {} HDUs'.format(image_filename.strip(), len(img_hdulist)))
        ccd_count += 1
    img_hdulist.close()

print(ccd_count, ccd_count/len(image_filename_unique))
