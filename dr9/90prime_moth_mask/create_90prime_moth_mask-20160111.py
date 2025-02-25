from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits
import healpy as hp
from astropy import wcs

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-90prime-dr9.fits.gz'
ccd = Table(fitsio.read(surveyccd_path))

# Only keep unique exposures
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]

expnum_min = 73990064
expnum_max = 73990184

mask = (ccd['expnum']>=expnum_min) & (ccd['expnum']<=expnum_max)
print(np.sum(mask))
expnum_list = np.unique(ccd['expnum'][mask])
print(len(expnum_list))
print(expnum_list)

for expnum in expnum_list:
    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)

    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3200:, 3100:] = True
    data['CCD1'] = mask

    np.savez_compressed(output_path, **data)

