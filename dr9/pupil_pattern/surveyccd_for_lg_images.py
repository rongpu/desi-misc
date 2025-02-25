from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

###########################################################################################################

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'

ccd_columns = ['image_filename']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
ccd['index'] = np.arange(len(ccd))

fns = [fn[fn.find('c4d_'):fn.find('_ooi_')+5] for fn in ccd['image_filename']]
fns = np.array(fns)

image_dir = '/global/cfs/cdirs/cosmo/staging/decam/CP-LG9'
image_path_list = sorted(glob.glob(os.path.join(image_dir, '*ooi*.fits.fz')))
print(len(image_path_list))
image_fn_list = [os.path.basename(image_path) for image_path in image_path_list]
image_fn_list1 = [image_fn[:image_fn.find('_ooi_')+5] for image_fn in image_fn_list]
image_fn_list1 = np.array(image_fn_list1)

mask = np.in1d(fns, image_fn_list1)
print(np.sum(mask))
ccd = ccd[mask]
fns = fns[mask]

ccd['image_filename_lg']=' '*len(image_fn_list[0])
for fn, fn1 in zip(image_fn_list, image_fn_list1):
    mask = fn1==fns
    ccd['image_filename_lg'][mask] = fn
    
ccd1 = Table(fitsio.read(surveyccd_path, rows=ccd['index']))
ccd1['image_filename_lg'] = ccd['image_filename_lg']

ccd1.write('/global/cfs/cdirs/desi/users/rongpu/dr9/pupil_pattern/survey-ccds-decam-dr9-lg.fits')

###########################################################################################################

blacklist = ['c4d_170722_103109_ooi_z_lg9.fits.fz', 'c4d_140813_094752_ooi_z_lg9.fits.fz', 'c4d_140814_095430_ooi_z_lg9.fits.fz', 'c4d_160226_063930_ooi_z_lg9.fits.fz', 'c4d_160227_064302_ooi_z_lg9.fits.fz', 'c4d_181002_051122_ooi_r_lg9.fits.fz', 'c4d_181002_051259_ooi_r_lg9.fits.fz', 'c4d_181101_033742_ooi_g_lg9.fits.fz', 'c4d_181102_034827_ooi_g_lg9.fits.fz', 'c4d_181102_035600_ooi_r_lg9.fits.fz', 'c4d_181103_030508_ooi_g_lg9.fits.fz', 'c4d_181104_035326_ooi_r_lg9.fits.fz', 'c4d_181104_035536_ooi_g_lg9.fits.fz', 'c4d_181104_040548_ooi_r_lg9.fits.fz', 'c4d_181104_045813_ooi_z_lg9.fits.fz', 'c4d_181104_050059_ooi_z_lg9.fits.fz', 'c4d_181105_043756_ooi_z_lg9.fits.fz']
ccd = Table.read('/global/cfs/cdirs/desi/users/rongpu/dr9/pupil_pattern/survey-ccds-decam-dr9-lg-gradient_added.fits')

ccd['blacklist'] = False
for image_filename_lg in blacklist:
    mask = ccd['image_filename_lg']==image_filename_lg
    ccd['blacklist'][mask] = True

ccd['run'] = -1
for band, run in zip(['g', 'r'], [1, 2]):
    mask = ccd['filter']==band
    ccd['run'][mask] = run

band = 'z'
mask = (ccd['filter']==band) & (ccd['mjd_obs']<58158.3519989+0.5) # cassette 3
ccd['run'][mask] = 3
mask = (ccd['filter']==band) & (ccd['mjd_obs']>=58158.3519989+0.5) # cassette 1
ccd['run'][mask] = 4

ccd.write('/global/cfs/cdirs/desi/users/rongpu/dr9/pupil_pattern/survey-ccds-decam-dr9-lg-pupilrun.fits', overwrite=True)

