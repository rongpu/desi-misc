from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

input_dir = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/calib/decam/psfex-merged'
# output_dir = '/global/cscratch1/sd/rongpu/temp/decam-psfex-merged'
output_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/new-psfex/decam-psfex-merged'

field = 'south'
region_name = 'decals_ngc'

radius_lim1, radius_lim2 = 4.5, 5.5
radius_lim3, radius_lim4 = 6., 6.8

ccd = fitsio.read('/global/project/projectdirs/cosmo/work/legacysurvey/dr9/reorg/survey-ccds-decam-dr9-newlocs2.fits.gz')
ccd = Table(ccd)
print(len(ccd))

# Load Schlegel's CCD file list
fn = '/global/project/projectdirs/cosmo/work/users/djschleg/dr9lists/dr9c.txt'
with open(fn, 'r') as f:
    lines = list(map(str.rstrip, f.readlines()))
print(len(lines))
# print(lines[0])

ccd['basename'] = list(map(os.path.basename, ccd['image_filename']))
mask = np.in1d(ccd['basename'], np.array(lines))
print(np.sum(mask)/len(mask))
ccd = ccd[mask]
print(len(ccd))

for band in ['g', 'r', 'z']:

    print('\n################################################################################################################')
    print('band = {}'.format(band))

    if (field=='north') and ((band=='g') or (band=='r')):
        pixscale_native = 0.454
    else:
        pixscale_native = 0.262
    # pixscale = 0.262 # pixscale for cutout queries

    ccd_mask = ccd['filter']==band
    print(np.sum(ccd_mask))

    expnum_all = np.sort(np.unique(ccd[ccd_mask]['expnum']))
    print('number of exposures:', len(expnum_all))

    tmp = np.loadtxt('/global/homes/r/rongpu/desi/star_profiles/{}_poly_fit.txt'.format(region_name))
    band_index = np.where(band==np.array(['g', 'r', 'z']))[0][0]
    poly = np.poly1d(tmp[band_index])
    print(poly)
    profile_fit = np.poly1d(poly)

    ################################################################################################################

    for index, expnum in enumerate(expnum_all):

        if (index+1)%(len(expnum_all)//10)==0:
            print('{:.0f}%'.format(index/(len(expnum_all))*100))

        expnum_str = str(expnum)

        input_fn = os.path.join(input_dir, '{}/decam-{}.fits'.format((5-len(expnum_str[:3]))*'0'+expnum_str[:3], (8-len(expnum_str))*'0'+expnum_str))
        output_fn = os.path.join(output_dir, '{}/decam-{}.fits'.format((5-len(expnum_str[:3]))*'0'+expnum_str[:3], (8-len(expnum_str))*'0'+expnum_str))
        
        # skip if the new psfex file already exists
        if os.path.isfile(output_fn):
            continue
        
        hdu = fits.open(input_fn)
        psf_mask = hdu[1].data['psf_mask']
        # print(psf_mask.shape)

        for ccd_index in range(len(psf_mask)):
            
            psf_all = psf_mask[ccd_index]
            # print(psf_all.shape)

            # psf_index = 0
            for psf_index in range(len(psf_all)):
                
                psfi = psf_all[psf_index]

                if psf_index==0:
                    # normalize to a 22.5 magnitude star
                    # print(np.sum(psfi))
                    psfi = psfi/np.sum(psfi)

                    grid = pixscale_native * np.linspace(-0.5*(psfi.shape[0]-1), 0.5*(psfi.shape[0]-1), psfi.shape[0])
                    xx, yy = np.meshgrid(grid, grid)
                    radius_grid = np.sqrt(xx**2 + yy**2)
                    radius = radius_grid.flatten()

                    with warnings.catch_warnings():
                        warnings.simplefilter("ignore")
                        img_fit = 10**(profile_fit(np.log10(radius_grid)))
                        img_fit[~np.isfinite(img_fit)] = 0

                    psfi_combine = psfi.copy()

                    r1, r2 = radius_lim1, radius_lim2
                    mask = (radius_grid>r1) & (radius_grid<r2)
                    psfi_combine[mask] = psfi[mask] * (r2-radius_grid[mask])/(r2-r1) \
                                   + img_fit[mask] * (radius_grid[mask]-r1)/(r2-r1)

                    r1, r2 = radius_lim2, radius_lim3
                    mask = (radius_grid>=r1) & (radius_grid<r2)
                    psfi_combine[mask] = img_fit[mask]

                    r1, r2 = radius_lim3, radius_lim4
                    mask = (radius_grid>=r1) & (radius_grid<r2)
                    psfi_combine[mask] = img_fit[mask] * (r2-radius_grid[mask])/(r2-r1) \
                                   + 0 * (radius_grid[mask]-r1)/(r2-r1)

                    mask = (radius_grid>radius_lim4)
                    psfi_combine[mask] = 0
                
                else:
                    
                    psfi_combine = psfi.copy()

                    r1, r2 = radius_lim1, radius_lim2
                    mask = (radius_grid>r1) & (radius_grid<r2)
                    psfi_combine[mask] = psfi[mask] * (r2-radius_grid[mask])/(r2-r1) \
                                   + 0 * (radius_grid[mask]-r1)/(r2-r1)
                    
                    mask = (radius_grid>radius_lim2)
                    psfi_combine[mask] = 0

                psf_all[psf_index] = psfi_combine

            psf_mask[ccd_index] = psf_all

        hdu[1].data['psf_mask'] = psf_mask

        if not os.path.exists(os.path.dirname(output_fn)):
            os.makedirs(os.path.dirname(output_fn))
        # hdu.writeto(output_fn)
        hdu.writeto(output_fn, overwrite=False)
