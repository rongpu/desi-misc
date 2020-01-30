# Use the pre-computed nightly sky model for the computation

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

from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'} 
plt.rcParams.update(params)


nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

ccdnamenumdict = {'S1': 25, 'S2': 26, 'S3': 27, 'S4':28,
                  'S5': 29, 'S6': 30, 'S7': 31,
                  'S8': 19, 'S9': 20, 'S10': 21, 'S11': 22, 'S12': 23,
                  'S13': 24,
                  'S14': 13, 'S15': 14, 'S16': 15, 'S17': 16, 'S18': 17,
                  'S19': 18,
                  'S20': 8, 'S21': 9, 'S22': 10, 'S23': 11, 'S24': 12,
                  'S25': 4, 'S26': 5, 'S27': 6, 'S28': 7,
                  'S29': 1, 'S30': 2, 'S31': 3,
                  'N1': 32, 'N2': 33, 'N3': 34, 'N4': 35,
                  'N5': 36, 'N6': 37, 'N7': 38,
                  'N8': 39, 'N9': 40, 'N10': 41, 'N11': 42, 'N12': 43,
                  'N13': 44,
                  'N14': 45, 'N15': 46, 'N16': 47, 'N17': 48, 'N18': 49,
                  'N19': 50,
                  'N20': 51, 'N21': 52, 'N22': 53, 'N23': 54, 'N24': 55,
                  'N25': 56, 'N26': 57, 'N27': 58, 'N28': 59,
                  'N29': 60, 'N30': 61, 'N31': 62,
                  }

fringe_dir = '/global/homes/d/djschleg/cosmo/staging/decam/DECam_CP-Fringe'
image_dir = '/global/project/projectdirs/cosmo/staging/'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'
blob_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'

sky_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/smooth_sky'
plot_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/plots'
output_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/data'

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'mjd_obs', 'ra', 'dec', 'skyrms', 'ccdraoff', 'ccddecoff', 'ccd_cuts']
ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
# ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
mask = ccd['ccd_cuts']==0
mask &= ccd['filter']=='z' # include only z-band images
ccd = ccd[mask]
print(len(ccd), 'CCD')

# Find CCDs around some MJD
mask = (ccd['mjd_obs']>(57815-4)) & (ccd['mjd_obs']<(57815+4)) # DECaLS observing run starting Feb 28, 2017
mask |= ((ccd['mjd_obs']>(58359-2)) & (ccd['mjd_obs']<(58359+27))) # Starting Aug 28, 2018
mask |= ((ccd['mjd_obs']>(58423-2)) & (ccd['mjd_obs']<(58423+30))) # Two runs starting Oct 28, 2018
mask |= ((ccd['mjd_obs']>(57893-2)) & (ccd['mjd_obs']<(57893+30))) # Two runs starting May 18, 2017
ccd = ccd[mask]
print(len(ccd))

img_list_all = []

for hdu_index in range(1, 62):

    # skip S7
    if hdu_index==31:
        continue

    print('hdu_index =', hdu_index)

    mask = ccd['image_hdu']==hdu_index
    ccd1 = ccd.copy()
    ccd1 = ccd1[mask]
    # print(len(ccd1))

    # Identify the observing date of each CCD
    str_loc = np.char.find(np.array(ccd1['image_filename'], dtype='str'), '/CP201')
    ccd1['obs_date'] = np.array([ccd1['image_filename'][i][str_loc[i]+1:str_loc[i]+11] for i in range(len(ccd1))])
    t = Table()
    t['date'], t['counts'] = np.unique(ccd1['obs_date'], return_counts=True)

    # Require a minimum number of CCDs (since scipy gaussian_filter does not handle NAN)
    mask = t['counts']<50
    mask_remove = np.in1d(ccd1['obs_date'], t['date'][mask])
    ccd1 = ccd1[~mask_remove]
    # print(len(ccd1))

    # Find the CCDs whose blobmask files exist
    ccd_mask = np.zeros(len(ccd1), dtype=bool)
    for ccd_index in range(len(ccd1)):
        str_loc = str.find(ccd1['image_filename'][ccd_index].strip(), '.fits')
        img_filename_base = ccd1['image_filename'][ccd_index].strip()[:str_loc]
        blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
        if os.path.isfile(blob_path):
            ccd_mask[ccd_index] = True
    # print(np.sum(ccd_mask)/len(ccd_mask))
    ccd1 = ccd1[ccd_mask]
    # print(len(ccd1))

    obs_date_list = np.unique(ccd1['obs_date'])
    print('Total nubmer of nights: ', len(obs_date_list))

    if os.path.isfile(os.path.join(output_dir, 'fringe_{}.npy'.format(hdu_index))):
        continue

    for obs_index, obs_date in enumerate(obs_date_list):

        ccd_mask = ccd1['obs_date']==obs_date
        print(obs_index, obs_date+':', np.sum(ccd_mask), 'CCDs')

        try:
            smooth_sky = np.load(os.path.join(sky_dir, 'smooth_sky_{}_{}.npy'.format(obs_date, hdu_index)))
        except:
            print(os.path.join(sky_dir, 'smooth_sky_{}_{}.npy'.format(obs_date, hdu_index)), 'does not exist!!!')
            continue

        for index, ccd_index in enumerate(np.where(ccd_mask)[0]):
            print(index)
            
            # Load CCD image
            img_fn = os.path.join(image_dir, ccd1['image_filename'][ccd_index]).strip()
            hdulist = fits.open(img_fn)
            # Some images do not have FRGSCALE in the header
            try:
                frgscale = (hdulist[ccd1['image_hdu'][ccd_index]].header)['FRGSCALE']
            except:
                print('frgscale does not exist!!!')
                continue
            # w = wcs.WCS(hdulist[ccd1['image_hdu'][ccd_index]].header)
            img = hdulist[ccd1['image_hdu'][ccd_index]].data
            
            # Load fringe image
            ccdnum = str(ccdnamenumdict[ccd1[ccd_index]['ccdname'].strip()]).zfill(2)
            fringe_path = os.path.join(fringe_dir, 'DES17B_20180103_908c062-z-{}_frg.fits'.format(ccdnum))
            fringe = fits.getdata(fringe_path)
            # remove the edge pixels
            fringe = fringe[1:4095, 1:2047]
            
            # Back out the exisiting fringe correction
            img += fringe*frgscale
            
            # Load blob mask
            str_loc = str.find(ccd1['image_filename'][ccd_index].strip(), '.fits')
            img_filename_base = ccd1['image_filename'][ccd_index].strip()[:str_loc]
            blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
            blob_data = np.load(blob_path)
            blob = blob_data['hdu'+str(hdu_index).zfill(2)]
            
            # Remove median sky
            sky = np.median(img[blob].flatten())
            img = img - sky
            
            # Normalize by frgscale
            img = img/frgscale
            
            # Apply blob mask
            img[~blob] = np.nan

            # # 3-sigma clipping
            # sky_nmad = nmad(img[np.isfinite(img)]) # sky level
            # mask = (img<-3*sky_nmad) | (img>3*sky_nmad)
            # img = img.copy()
            # img[mask] = np.nan

            # Subtract off the smooth component
            img -= smooth_sky

            img_list_all.append(img)

        print('=============================================================================')

    vrange = 5e-3

    img_median_final = np.nanmedian(img_list_all, axis=0)
    np.save(os.path.join(output_dir, 'fringe_{}.npy'.format(hdu_index)), img_median_final)

    plt.figure(figsize=(17, 8))
    plt.imshow((img_median_final).T, cmap='seismic', vmin=-vrange, vmax=vrange)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'final_fringe_{}.png'.format(hdu_index)))
    plt.close()

    # Plot 4-pixel gaussian smoothed fringe image
    # 5-sigma clipping
    sky_nmad = nmad(img_median_final[np.isfinite(img_median_final)]) # sky level
    img_median_final1 = img_median_final.copy()
    # mask = (img_median_final<-5*sky_nmad) | (img_median_final>5*sky_nmad)
    # img_median_final1[mask] = 0
    mask = (img_median_final<-5*sky_nmad)
    img_median_final1[mask] = -5*sky_nmad
    mask = (img_median_final>5*sky_nmad)
    img_median_final1[mask] = 5*sky_nmad
    img_median_4pix_gauss = gaussian_filter((img_median_final1), 4, mode='reflect')
    plt.figure(figsize=(17, 8))
    plt.imshow((img_median_4pix_gauss).T, cmap='seismic', vmin=-vrange, vmax=vrange)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'final_fringe_smooth_{}.png'.format(hdu_index)))
    plt.close()

    plt.figure(figsize=(17, 8))
    plt.imshow((fringe).T, cmap='seismic', vmin=-vrange, vmax=vrange)
    plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, 'original_fringe_{}.png'.format(hdu_index)))
    plt.close()

