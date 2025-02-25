from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

ccdnum_list = [1, 2, 3, 4]
bit_to_set = 1

image_dir = '/global/cfs/cdirs/cosmo/staging/'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/90prime_moth_mask/20160603'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-90prime-dr9.fits.gz'

ccd = Table(fitsio.read(surveyccd_path))

mask = ccd['ccd_cuts']==0
ccd = ccd[mask]

# Only keep unique exposures
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]

expnum_list = [75430048, 75430049, 75430050, 75430051, 75430052, 75430053, 75430054, 75430055, 75430056, 75430057, 75430073, 75430074, 75430075, 75430076, 75430077, 75430078, 75430079, 75430080, 75430082, 75430083, 75430084, 75430087, 75430088, 75430089, 75430090, 75430091, 75430092, 75430093, 75430094, 75430095, 75430096, 75430097, 75430098, 75430099, 75430100, 75430101, 75430102, 75430103, 75430104, 75430105, 75430106, 75430107, 75430108, 75430109, 75430110, 75430111, 75430112, 75430113, 75430114, 75430115, 75430116, 75430117, 75430118, 75430119, 75430120, 75430121, 75430122, 75430123, 75430124, 75430125, 75430126]
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
