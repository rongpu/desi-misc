from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

ccdnum_list = [1, 2, 3, 4]
bit_to_set = 1

image_dir = '/global/cfs/cdirs/cosmo/staging/'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/90prime_moth_mask/20170619'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-90prime-dr9.fits.gz'

ccd = Table(fitsio.read(surveyccd_path))

mask = ccd['ccd_cuts']==0
ccd = ccd[mask]

# Only keep unique exposures
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]

expnum_list = [79240047, 79240048, 79240049, 79240050, 79240057, 79240058, 79240059, 79240060, 79240061, 79240062, 79240063, 79240064, 79240065]
print(len(expnum_list))

idx = np.where(np.in1d(ccd['expnum'], expnum_list))[0]
if len(idx)!=len(expnum_list):
    raise ValueError

for ccd_index in idx:

    fn = ccd['image_filename'][ccd_index].strip().replace('_ooi_', '_ood_')
    print(fn)
    image_path = os.path.join(image_dir, fn)
    output_path = os.path.join(output_dir, fn)
    if not os.path.isdir(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))

    expnum = ccd['expnum'][ccd_index]
    band = ccd['filter'][ccd_index]
    print(ccd_index, band, expnum)

    mask_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = np.load(mask_path)
    hdul = fitsio.FITS(output_path, mode='rw', clobber=True)
    hdr = fitsio.read_header(image_path, ext=0)
    hdul.write(data=None, header=hdr)  # first HDU is empty

    for ii, ccdnum in enumerate(ccdnum_list):

        ccdname = 'CCD'+str(ccdnum)
        ccdname_lowercase = 'ccd'+str(ccdnum)
        img = fitsio.read(image_path, ext=ccdname_lowercase)
        hdr = fitsio.read_header(image_path, ext=ccdname_lowercase)

        if ccdname in data.keys():
            new_mask = data[ccdname]
            mask = new_mask & (img & (2**bit_to_set)==0)
            img[mask] += 2**bit_to_set

        hdul.write(data=img, header=hdr, extname=ccdname_lowercase, compress='rice')

    hdul.close()
