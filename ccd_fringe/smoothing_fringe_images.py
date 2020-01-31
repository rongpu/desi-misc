from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits
import healpy as hp
from astropy import wcs

from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'} 
plt.rcParams.update(params)

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

fringe_dir = '/global/homes/d/djschleg/cosmo/staging/decam/DECam_CP-Fringe'
image_dir = '/global/project/projectdirs/cosmo/staging/'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'
blob_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'

sky_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/smooth_sky'
plot_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/plots_postproc'
data_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/data'

# CP mask image from an arbitrary exposure
ood_fn = '/global/project/projectdirs/cosmo/staging/decam/CP/V4.8.2/CP20170327/c4d_170328_100215_ood_z_ls9.fits.fz'


for ii, hdu_index in enumerate(range(1, 62)):
    
    if hdu_index==31:
        print('skipping S7')
        # sys.exit()
        continue

    print(hdu_index)

    fringe = np.load(os.path.join(data_dir, 'fringe_{}.npy'.format(hdu_index)))
    hdulist = fits.open(ood_fn)
    ood = hdulist[hdu_index].data
    print('Fraction of masked pixels:', np.sum(ood!=0)/len(ood.flatten()))

    fringe_masked = fringe.copy()
    fringe_masked[ood!=0] = np.nan

    sky_nmad = nmad(fringe_masked[np.isfinite(fringe_masked)])

    mask = (fringe_masked<-3*sky_nmad) | (fringe_masked>3*sky_nmad)
    print('Fraction of clipped pixels:', np.sum(mask)/np.sum(np.isfinite(fringe_masked)))
    fringe_masked[mask] = np.nan
    
    vrange = 5e-3

    plt.figure(figsize=(17, 8))
    plt.imshow((fringe_masked).T, cmap='seismic', vmin=-vrange, vmax=vrange)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'new_fringe_masked_{}.png'.format(hdu_index)))
    plt.close()

    kernel = Gaussian2DKernel(stddev=2)
    fringe_smooth = convolve(fringe_masked, kernel, boundary='extend', nan_treatment='interpolate')

    plt.figure(figsize=(17, 8))
    plt.imshow((fringe_smooth).T, cmap='seismic', vmin=-vrange, vmax=vrange)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'new_fringe_smooth_{}.png'.format(hdu_index)))
    plt.close()

    np.save(os.path.join(data_dir, 'fringe_smooth_{}.npy'.format(hdu_index)), fringe_smooth)

