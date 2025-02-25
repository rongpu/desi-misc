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

from multiprocessing import Pool

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

# plot_dir = '/global/cscratch1/sd/rongpu/fringe/plots'
# output_dir = '/global/cscratch1/sd/rongpu/fringe/smooth_sky'

plot_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/plots'
output_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/smooth_sky'

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'ccdname', 'filter', 'mjd_obs', 'ccd_cuts']
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

def check_progress(hdu_index):

    # skip S7
    if hdu_index==31:
        pass

    # print(hdu_index)

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

    ##############################################################################################################

    # # Compute the median stacked image for one specific night
    # np.random.seed(123)
    # obs_date = np.random.choice(ccd1['obs_date'])
    # ccd_mask = ccd1['obs_date']==obs_date
    # print(obs_date+':', np.sum(ccd_mask), 'CCDs')

    obs_date_list, n_exp = np.unique(ccd1['obs_date'], return_counts=True)
    obs_date_list = obs_date_list[np.argsort(n_exp)]
    # print('Total nubmer of nights: ', len(obs_date_list))

    # # Randomly select a few nights
    # obs_date_list = np.array(obs_date_list)
    # np.random.seed(123)
    # obs_date_list = np.random.choice(obs_date_list, size=10, replace=False)

    obs_date_complete = np.zeros(len(obs_date_list), dtype=bool)

    for obs_index, obs_date in enumerate(obs_date_list):
        if os.path.isfile(os.path.join(output_dir, 'smooth_sky_{}_{}.npy'.format(obs_date, hdu_index))):
            obs_date_complete[obs_index] = True

    print(hdu_index, len(obs_date_list), '{:.0f}%'.format(np.sum(obs_date_complete)/len(obs_date_complete)*100))

for hdu_index in np.arange(1, 62):
    check_progress(hdu_index)

