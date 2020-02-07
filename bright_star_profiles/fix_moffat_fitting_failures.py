from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'
psfex_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/calib/patched-psfex'
psfex_dir_new = '/global/project/projectdirs/desi/users/rongpu/dr9/patched-psfex'

ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts']
ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
ccd = Table(ccd)
print(len(ccd))

mask = ccd['ccd_cuts']==0
print(np.sum(mask)/len(mask))
ccd = ccd[mask]

unique_expnum = np.unique(ccd['expnum'])
print(len(unique_expnum))

exp_index_list = np.arange(len(unique_expnum))

exp_index_with_zero_values = []
exp_index_no_exist = []

for exp_index in exp_index_list:

    if exp_index%(len(exp_index_list)//100)==0:
        print(exp_index, '/', len(exp_index_list))

    mask = ccd['expnum']==unique_expnum[exp_index]
    band = ccd['filter'][mask][0]
    # print('band = {}'.format(band))

    image_filename = ccd['image_filename'][mask][0]
    psfex_filename = image_filename[:image_filename.find('.fits.fz')]+'-psfex.fits'
    psfex_filename_new = image_filename[:image_filename.find('.fits.fz')]+'-psfex-patched.fits'
    psfex_path = os.path.join(psfex_dir, psfex_filename)
    
    # print(psfex_path)
    
    if not os.path.isfile(psfex_path):
        exp_index_no_exist.append(exp_index)
        print(exp_index, 'does not exist')
        continue

    data = fitsio.read(psfex_path, columns=['moffat_alpha', 'moffat_beta'])
    data = Table(data)
        
    if (not np.all(data['moffat_alpha']!=0)) or (not np.all(data['moffat_beta']!=0)):
        
        exp_index_with_zero_values.append(exp_index)
        print(exp_index, band+'-band invalid moffat parameters')

        if not np.all(data['moffat_alpha']==0):
            mask = data['moffat_alpha']!=0
            moffat_alpha_median = np.median(data['moffat_alpha'][mask])
            moffat_beta_median = np.median(data['moffat_beta'][mask])
            
            data['moffat_alpha'][~mask] = moffat_alpha_median
            data['moffat_beta'][~mask] = moffat_beta_median
            print(exp_index, 'median_alpha = {:.4f},  median_beta = {:.4f}'.format(moffat_alpha_median, moffat_beta_median))
        else:
            print(exp_index, 'using default values')
            data['moffat_alpha'] = 0.8
            data['moffat_beta'] = 2.2

        psfex_path_new = os.path.join(psfex_dir_new, psfex_filename)
        if not os.path.exists(os.path.dirname(psfex_path_new)):
            os.makedirs(os.path.dirname(psfex_path_new))
        data.write(psfex_path_new)
