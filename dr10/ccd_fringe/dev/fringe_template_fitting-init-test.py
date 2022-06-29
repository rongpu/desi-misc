# Compute the initial fringe templates to get the initial fringe scales for normalization

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio

from pathlib import Path
from multiprocessing import Pool

from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter
from scipy import stats


time_start = time.time()

params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

nmad = lambda x: 1.4826 * np.nanmedian(np.abs(x-np.nanmedian(x)))

ccdnamenumdict = {'S1': 25, 'S2': 26, 'S3': 27, 'S4': 28,
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
ccdnamenumdict_inv = {aa: bb for bb, aa in ccdnamenumdict.items()}

img_shape = (4094, 2046)

# Trim CCD edges
img_trim_size_x = 15 + 224
img_trim_size_y = 15 + 224

fringe_dr9_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/fringe/DECam_CP-Fringe-Normed'
fringe_init_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/fringe_templates_init'
image_dir = '/global/cfs/cdirs/cosmo/staging/'
surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/survey-ccds-z-fringe.fits'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_ccd_blob_mask'

output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/frgscale_init_test'

# Load DR9 fringe templates
fringe_templates_dr9 = {}
for ccdnum in range(1, 63):
    ccdname = ccdnamenumdict_inv[ccdnum]
    fringe_template_path = os.path.join(fringe_dr9_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
    if os.path.isfile(fringe_template_path):
        fringe_tmp = fitsio.read(fringe_template_path)
        fringe_tmp = fringe_tmp[1:4095, 1:2047]  # remove the edge pixels
        fringe_tmp = fringe_tmp[img_trim_size_y:(fringe_tmp.shape[0]-img_trim_size_y), img_trim_size_x:(fringe_tmp.shape[1]-img_trim_size_x)]
        fringe_templates_dr9[ccdnum] = fringe_tmp.copy()

# Load initial fringe templates
fringe_templates_init = {}
for ccdnum in range(1, 63):
    ccdname = ccdnamenumdict_inv[ccdnum]
    fringe_template_path = os.path.join(fringe_init_dir, 'DECam_z_frg_{}_CCD{}_smooth.fits'.format(ccdname, str(ccdnum).zfill(2)))
    if os.path.isfile(fringe_template_path):
        fringe_tmp = fitsio.read(fringe_template_path)
        fringe_tmp = fringe_tmp[1:4095, 1:2047]  # remove the edge pixels
        fringe_tmp = fringe_tmp[img_trim_size_y:(fringe_tmp.shape[0]-img_trim_size_y), img_trim_size_x:(fringe_tmp.shape[1]-img_trim_size_x)]
        fringe_templates_init[ccdnum] = fringe_tmp.copy()

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'ccdname', 'filter', 'expnum', 'exptime', 'ccdskycounts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
print(len(ccd))

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx].copy()
str_loc = np.char.find(np.array(exp['image_filename'], dtype='str'), '/CP20')
exp['obs_date'] = np.array([exp['image_filename'][i][str_loc[i]+3:str_loc[i]+11] for i in range(len(exp))])
print('Total number of nights:', len(np.unique(exp['obs_date'])))
print('Number of images: {}'.format(len(exp)))

ccd_all = ccd.copy()

overwrite = True


def fringe_template_fitting(expnum):

    print('expnum:', expnum)
    mask = exp['expnum']==expnum
    image_fn = exp['image_filename'][mask][0]

    mask = ccd_all['expnum']==expnum
    ccd = ccd_all[mask].copy()

    frgscale_output_path = os.path.join(output_dir, image_fn.strip().replace('.fits.fz', '.txt'))

    if os.path.isfile(frgscale_output_path) and not overwrite:
        print(frgscale_output_path+' already exists!')
        return None

    if not os.path.exists(os.path.dirname(frgscale_output_path)):
        try:
            os.makedirs(os.path.dirname(frgscale_output_path))
        except:
            pass

    fringe_params = []

    for ccd_index in range(len(ccd)):

        ccdname = ccd['ccdname'][ccd_index]
        ccdnum = ccdnamenumdict[ccdname]

        ################################################
        if not ccdnum in fringe_templates_init.keys():
            # print(ccdnum, 'not in fringe_templates_init')
            continue
        ################################################

        # skip S7 and S30
        if (ccdnum==2) or (ccdnum==31):
            continue

        hdu_index = ccd['image_hdu'][ccd_index]

        # Load CCD image
        img_fn = os.path.join(image_dir, image_fn).strip()
        ood_fn = img_fn.replace('_ooi_', '_ood_')

        header = fitsio.read_header(img_fn, ext=hdu_index)
        if ('FRGSCNEW' not in header.keys()) and ('FRGSCALE' in header.keys()):
            raise ValueError(expnum, ccdname, image_fn)
        if 'FRGSCNEW' in header.keys():
            frgscale = header['FRGSCNEW']
        else:
            # print(expnum, ccdname, 'FRGSCNEW does not exist!!!')
            frgscale = 0
        # w = wcs.WCS(hdulist[ccd['image_hdu'][ccd_index]].header)

        img = fitsio.read(img_fn, ext=hdu_index)
        ood = fitsio.read(ood_fn, ext=hdu_index)

        # Back out the exisiting fringe correction
        if frgscale!=0:
            img += fringe_templates_dr9[ccdnum] * frgscale

        # Load blob mask
        str_loc = str.find(image_fn.strip(), '.fits')
        img_filename_base = image_fn.strip()[:str_loc]
        blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
        try:
            blob_data = np.load(blob_path)
        except:
            print('Error loading {}'.format(blob_path))
            continue
        try:
            blob = blob_data['hdu'+str(hdu_index).zfill(2)]
        except:
            print('Error loading {} blobmask {}'.format(ccdname, blob_path))
            continue

        mask_good = blob & (ood==0)

        # Remove median sky
        median_sky = np.median(img[mask_good].flatten())
        img = img - median_sky

        # Apply mask
        img[~mask_good] = np.nan

        img = img[img_trim_size_y:(img.shape[0]-img_trim_size_y), img_trim_size_x:(img.shape[1]-img_trim_size_x)]
        print(img.shape)

        if np.all(np.isnan(img)):
            print('{} CCD {} are all NaNs'.format(img_fn, ccdname))

        # 3-sigma clipping
        sky_nmad = nmad(img)
        mask = (img<-3*sky_nmad)
        img[mask] = -3*sky_nmad
        mask = (img>3*sky_nmad)
        img[mask] = 3*sky_nmad

        #################################### Subtract smooth sky ####################################

        img1 = img.copy()

        # # trim edges
        # trim_size = 15   # trimmed imagee size need to be multiples of binsize
        # img1 = img1[trim_size:(img1.shape[0]-trim_size), trim_size:(img1.shape[1]-trim_size)]
        trim_size = 0

        # downsize the image to speed up gaussian filter
        binsize = 32
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            img1 = np.nanmedian(np.nanmedian(img1.reshape((img1.shape[0]//binsize, binsize, img1.shape[1]//binsize,-1)), axis=3), axis=1)
        x_small_grid = trim_size + binsize/2+binsize*np.arange(img1.shape[1])
        y_small_grid = trim_size + binsize/2+binsize*np.arange(img1.shape[0])
        img1[np.isnan(img1)] = 0

        # Gaussian filtering
        img1_smooth = gaussian_filter(img1, 3, mode='reflect')

        # Convert the downsized smooth image to full size
        interp_func = interp2d(x_small_grid, y_small_grid, img1_smooth, kind='linear')
        x_grid, y_grid = np.arange(img.shape[1]), np.arange(img.shape[0])
        smooth_sky = interp_func(x_grid, y_grid).reshape(img.shape)

        # Subtract off the smooth component
        img -= smooth_sky

        #############################################################################################

        # Linear regression
        img_mask = np.isfinite(img)
        img_mask &= np.isfinite(fringe_templates_init[ccdnum])
        if np.sum(img_mask)==0:
            continue
        img_flat = img[img_mask].flatten()
        fringe_flat = fringe_templates_init[ccdnum][img_mask].flatten()
        slope, intercept, r_value, p_value, std_err = stats.linregress(fringe_flat, img_flat)

        fringe_params.append([expnum, ccdnum, ccdname, frgscale, median_sky, sky_nmad, slope, intercept, r_value, p_value, std_err])

    if len(fringe_params)==0:
        Path(frgscale_output_path).touch()
    else:
        fringe_table = Table(np.array(fringe_params), names=['expnum', 'ccdnum', 'ccdname', 'frgscale_old', 'median_sky', 'sky_nmad', 'slope', 'intercept', 'r_value', 'p_value', 'std_err'])
        fringe_table['expnum'] = np.array(fringe_table['expnum'], dtype=int)
        fringe_table['ccdnum'] = np.array(fringe_table['ccdnum'], dtype=int)
        fringe_table.write(frgscale_output_path, format='ascii.commented_header', overwrite=True)

    return None


print('=============================================================================')

# expnum_list = np.array(exp['expnum'])

# n_processes = 32
# with Pool(processes=n_processes) as pool:
#     res = pool.map(fringe_template_fitting, expnum_list, chunksize=1)

expnum = 865626
time_start = time.time()
fringe_template_fitting(expnum)
print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))

