# Create high-resolution images.

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from scipy.ndimage.filters import gaussian_filter

ccdnum_list = [1, 2, 3, 4]
ccd_ra = [0.2813, 0.2813, -0.2813, -0.2813]
ccd_dec = [0.263, -0.263, 0.263, -0.263]

# surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-90prime-dr9.fits.gz'
# surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-mosaic-dr9.fits.gz'

ccd = Table(fitsio.read(surveyccd_path))

# # Only keep unique exposures
# _, idx = np.unique(ccd['expnum'], return_index=True)
# ccd = ccd[idx]

###############################################################

expnum = 82310056
ccd_index = np.where(ccd['expnum']==expnum)[0][0]

binsize = 4
pix_size = 0.454/3600*binsize
vrange = 1.

# expnum = ccd['expnum'][ccd_index]
band = ccd['filter'][ccd_index]
print(ccd_index, band, expnum)
fn = ccd['image_filename'][ccd_index]

plt.figure(figsize=(14, 14))

for ii, ccdnum in enumerate(ccdnum_list):
    
    try:
        img = fits.getdata('/global/cfs/cdirs/cosmo/staging/'+fn, extname='CCD'+str(ccdnum))
    except:
        print('Failure loading {}'.format('/global/cfs/cdirs/cosmo/staging/'+fn))
        continue

    ################ downsize image ################

    if binsize!=1:

        trim_size_x = img.shape[1] % binsize
        trim_size_y = img.shape[0] % binsize
        img = img[:(img.shape[0]-trim_size_y), :(img.shape[1]-trim_size_x)]
        # to ignore NAN values, use np.nanmean
        img = np.nanmean(np.nanmean(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)

    ################################################

    ysize, xsize = img.shape
    ra, dec = ccd_ra[ii], ccd_dec[ii]

    # naive sky estimation
    mask = (img<np.percentile(img.flatten(), 95))
    median_sky = np.median(img[mask].flatten())
    img = img - median_sky

    img[~np.isfinite(img)] = 0
    # img = gaussian_filter(img, 2, mode='reflect', truncate=3)
    fig = plt.imshow(img.T, cmap='seismic', vmin=-vrange, vmax=vrange, 
               extent=(ra+ysize*pix_size/2, ra-ysize*pix_size/2, dec-xsize*pix_size/2, dec+xsize*pix_size/2))

plt.axis([0.55, -0.55, -0.55, 0.55])
plt.axis('off')
fig.axes.get_xaxis().set_visible(False)
fig.axes.get_yaxis().set_visible(False)
plt.colorbar(fraction=0.04, pad=0.04)
plt.tight_layout()
plt.show()

