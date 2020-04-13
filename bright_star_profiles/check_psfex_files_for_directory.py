from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

output_dir = '/global/cfs/projectdirs/cosmo/work/legacysurvey/dr9/calib/patched-psfex/decam/CP/V4.1'
# output_dir = '/global/cfs/projectdirs/cosmo/work/legacysurvey/dr9/calib/v2-patched-psfex'
# surveyccd_path = '/global/cfs/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'

file_list = glob.glob(os.path.join(output_dir, '*/*.fits'))
print(len(file_list))
print(file_list[0])

files_to_reprocess = []

for index in range(len(file_list)):

    if index%100==0:
        print(index, '/', len(file_list))

    psfex_path_new = file_list[index]

    with fitsio.FITS(psfex_path_new) as hdu:
        if not 'moffat_alpha' in hdu[1].get_colnames():
            files_to_reprocess.append(index)
            # print('found one: expnum =', expnum)
            print(psfex_path_new)

if len(files_to_reprocess)>0:
    pass
