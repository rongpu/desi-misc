from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

##########################################################################################

ccd = Table.read('/global/cscratch1/sd/rongpu/dr9dev/pupil_pattern/survey-ccds-decam-dr9-pupil-params.fits')
skyrun = Table.read('/global/cscratch1/sd/schlafly/legacysurvey/skyrunsgoodcountexpnumv48dr8.fits')
skyrun.sort('expnum')
ccd.sort('expnum')
if not np.sum(np.in1d(ccd['expnum'], skyrun['expnum']))==len(skyrun):
    raise ValueError
ccd['good_ccds'] = -1
mask = np.in1d(ccd['expnum'], skyrun['expnum'])
ccd['good_ccds'][mask] = skyrun['good_ccds']

mask = np.isfinite(ccd['PUPILAMP'])
print(np.sum(~mask))
ccd = ccd[mask]

mask = ccd['good_ccds']>10
ccd = ccd[mask]
print(len(ccd))

n_plots = 3
percentiles = [0, 25, 50, 75, 90, 95, 96, 97, 97.5, 98, 98.5, 99, 99.25, 99.5, 99.75, 100]

for band in ['g', 'r', 'z']:

    print(band)
    mask_band = ccd['filter']==band
    values = np.percentile(ccd['PUPILAMP'][mask_band], percentiles)
    print(values)
    
    f = open("/global/cfs/cdirs/desi/users/rongpu/plots/dr9dev/pupil_pattern/{}-band-images-CP.html".format(band), "w")
    f.write('<html>\n')
    f.write('<table>\n')

    for pupilamp in values:
        idx = np.argsort(np.abs(ccd['PUPILAMP'][mask_band]-pupilamp))
        expnum_list = list(ccd['expnum'][mask_band][idx[:3]])

        f.write('<tr>\n')
        for expnum in expnum_list:

            ccd_index = np.where(ccd['expnum']==expnum)[0][0]
            image_fn = os.path.basename(ccd['image_filename'][ccd_index].replace('.fits.fz', '.png'))

            f.write('<td><a href=\'cp_pupil/{}\'><img src=\'cp_pupil/{}\' width=\'400\'></a></td>\n'.format(image_fn, image_fn))
            # f.write('<td><a href=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\'><img src=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\' width=\'400\'></a></td>\n')
            
        f.write('</tr>\n')
            
    f.write('</table>\n')
    f.close()

##########################################################################################

ccd = Table.read('/global/cscratch1/sd/rongpu/dr9dev/pupil_pattern/survey-ccds-decam-dr9-pupil-params.fits')
skyrun = Table.read('/global/cscratch1/sd/schlafly/legacysurvey/skyrunsgoodcountexpnumv48dr8.fits')
skyrun.sort('expnum')
ccd.sort('expnum')
if not np.sum(np.in1d(ccd['expnum'], skyrun['expnum']))==len(skyrun):
    raise ValueError
ccd['good_ccds'] = -1
mask = np.in1d(ccd['expnum'], skyrun['expnum'])
ccd['good_ccds'][mask] = skyrun['good_ccds']

mask = np.isfinite(ccd['PUPILAMP'])
print(np.sum(~mask))
ccd = ccd[mask]

mask = ccd['good_ccds']>10
ccd = ccd[mask]
print(len(ccd))

n_plots = 3
percentiles = [0, 25, 50, 75, 90, 95, 96, 97, 97.5, 98, 98.5, 99, 99.25, 99.5, 99.75, 100]

for band in ['g', 'r', 'z']:

    print(band)
    mask_band = ccd['filter']==band
    values = np.percentile(ccd['PUPILMAX'][mask_band], percentiles)
    print(values)
    
    f = open("/global/cfs/cdirs/desi/users/rongpu/plots/dr9dev/pupil_pattern/{}-band-images-CP-pupilmax.html".format(band), "w")
    f.write('<html>\n')
    f.write('<table>\n')

    for pupilamp in values:
        idx = np.argsort(np.abs(ccd['PUPILMAX'][mask_band]-pupilamp))
        expnum_list = list(ccd['expnum'][mask_band][idx[:3]])

        f.write('<tr>\n')
        for expnum in expnum_list:

            ccd_index = np.where(ccd['expnum']==expnum)[0][0]
            image_fn = os.path.basename(ccd['image_filename'][ccd_index].replace('.fits.fz', '.png'))

            f.write('<td><a href=\'cp_pupilmax/{}\'><img src=\'cp_pupilmax/{}\' width=\'400\'></a></td>\n'.format(image_fn, image_fn))
            # f.write('<td><a href=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\'><img src=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\' width=\'400\'></a></td>\n')
            
        f.write('</tr>\n')
            
    f.write('</table>\n')
    f.close()
