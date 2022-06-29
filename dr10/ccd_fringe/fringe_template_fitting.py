# Get the fringe scales for the DECam deep exposures

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
fringe_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/fringe_templates'
image_dir = '/global/cfs/cdirs/cosmo/staging/'
# surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/survey-ccds-z-fringe.fits'
surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1.fits'

blob_dir_dr9 = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
blob_dir_dr10 = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_ccd_blob_mask'

output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/frgscale'

# Load DR9 fringe templates
fringe_templates_dr9 = {}
for ccdnum in range(1, 63):
    ccdname = ccdnamenumdict_inv[ccdnum]
    fringe_template_path = os.path.join(fringe_dr9_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
    if os.path.isfile(fringe_template_path):
        fringe_tmp = fitsio.read(fringe_template_path)
        fringe_tmp = fringe_tmp[1:4095, 1:2047]  # remove the edge pixels
        # fringe_tmp = fringe_tmp[img_trim_size_y:(fringe_tmp.shape[0]-img_trim_size_y), img_trim_size_x:(fringe_tmp.shape[1]-img_trim_size_x)]
        fringe_templates_dr9[ccdnum] = fringe_tmp.copy()

# Load initial fringe templates
fringe_templates = {}
for ccdnum in range(1, 63):
    ccdname = ccdnamenumdict_inv[ccdnum]
    fringe_template_path = os.path.join(fringe_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
    if os.path.isfile(fringe_template_path):
        fringe_tmp = fitsio.read(fringe_template_path)
        fringe_tmp = fringe_tmp[1:4095, 1:2047]  # remove the edge pixels
        fringe_tmp = fringe_tmp[img_trim_size_y:(fringe_tmp.shape[0]-img_trim_size_y), img_trim_size_x:(fringe_tmp.shape[1]-img_trim_size_x)]
        fringe_templates[ccdnum] = fringe_tmp.copy()
print('fringe_templates', len(fringe_templates))

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'ccdname', 'filter', 'expnum', 'exptime', 'ccdskycounts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
print(len(ccd))

mask = ccd['filter']=='z'
ccd = ccd[mask]
print(len(ccd))

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx].copy()
str_loc = np.char.find(np.array(exp['image_filename'], dtype='str'), '/CP20')
exp['obs_date'] = np.array([exp['image_filename'][i][str_loc[i]+3:str_loc[i]+11] for i in range(len(exp))])
print('Total number of nights:', len(np.unique(exp['obs_date'])))
print('Number of images: {}'.format(len(exp)))

overwrite = False


def fringe_template_fitting(expnum):

    print('expnum:', expnum)
    mask = exp['expnum']==expnum
    image_filename = exp['image_filename'][mask][0]
    img_fn = os.path.join(image_dir, image_filename).strip()
    ood_fn = img_fn.replace('_ooi_', '_ood_')

    frgscale_output_path = os.path.join(output_dir, image_filename.strip().replace('.fits.fz', '.txt'))

    n_hdu = len(fitsio.FITS(img_fn))

    if os.path.isfile(frgscale_output_path) and (not overwrite):
        tmp = Table.read(frgscale_output_path, format='ascii.commented_header')
        # if len(tmp)==n_hdu-1:
        print(frgscale_output_path+' already exists!')
        return None

    if not os.path.exists(os.path.dirname(frgscale_output_path)):
        try:
            os.makedirs(os.path.dirname(frgscale_output_path))
        except:
            pass

    fringe_params = []

    for hdu_index in range(1, n_hdu):

        header = fitsio.read_header(img_fn, ext=hdu_index)
        ccdname = header['EXTNAME'].strip()
        ccdnum = ccdnamenumdict[ccdname]

        # if not ccdnum in fringe_templates.keys():
        #     # print(ccdnum, 'not in fringe_templates')
        #     continue

        # # skip S7 and S30
        # if (ccdnum==2) or (ccdnum==31):
        #     continue

        if ('FRGSCNEW' not in header.keys()) and ('FRGSCALE' in header.keys()):
            raise ValueError(expnum, ccdname, img_fn)
        if 'FRGSCNEW' in header.keys():
            frgscale_dr9 = header['FRGSCNEW']
        else:
            # print(expnum, ccdname, 'FRGSCNEW does not exist!!!')
            frgscale_dr9 = 0
        # w = wcs.WCS(hdulist[ccd['image_hdu'][ccd_index]].header)

        img = fitsio.read(img_fn, ext=hdu_index)
        ood = fitsio.read(ood_fn, ext=hdu_index)

        # Back out the exisiting fringe correction
        if frgscale_dr9!=0:
            img += fringe_templates_dr9[ccdnum] * frgscale_dr9

        # Load blob mask
        str_loc = str.find(image_filename.strip(), '.fits')
        img_filename_base = image_filename.strip()[:str_loc]

        # Try DR9 blob mask first
        blob_path = os.path.join(blob_dir_dr9, 'blob_mask', img_filename_base+'-blobmask.npz')
        if (not os.path.isfile(blob_path)) or (os.stat(blob_path).st_size==0):
            blob_path = os.path.join(blob_dir_dr10, 'blob_mask', img_filename_base+'-blobmask.npz')

        if not os.path.isfile(blob_path):
            print('Blob mask does not exist', blob_path)
            continue
        elif (os.stat(blob_path).st_size==0):
            print('Blob mask is empty', blob_path)
            continue
        else:
            blob_data = np.load(blob_path)

        hdu_name = 'hdu'+str(hdu_index).zfill(2)
        if hdu_name in blob_data.files:
            blob = blob_data[hdu_name]
        else:
            print('{} ({}) does not exist in blobmask {}'.format(hdu_name, ccdname, blob_path))
            continue

        mask_good = blob & (ood==0)

        # Remove median sky
        median_sky = np.median(img[mask_good].flatten())
        img = img - median_sky

        # Apply mask
        img[~mask_good] = np.nan

        img = img[img_trim_size_y:(img.shape[0]-img_trim_size_y), img_trim_size_x:(img.shape[1]-img_trim_size_x)]

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
        img_mask &= np.isfinite(fringe_templates[ccdnum])
        if np.sum(img_mask)==0:
            continue
        img_flat = img[img_mask].flatten()
        fringe_flat = fringe_templates[ccdnum][img_mask].flatten()
        slope, intercept, r_value, p_value, std_err = stats.linregress(fringe_flat, img_flat)
        frgscale = slope

        fringe_params.append([expnum, ccdnum, ccdname, median_sky, sky_nmad, intercept, r_value, p_value, std_err, frgscale_dr9, frgscale])

    if len(fringe_params)==0:
        Path(frgscale_output_path).touch()
    else:
        fringe_table = Table(np.array(fringe_params), names=['expnum', 'ccdnum', 'ccdname', 'median_sky', 'sky_nmad', 'intercept', 'r_value', 'p_value', 'std_err', 'frgscale_dr9', 'frgscale'])
        fringe_table['expnum'] = np.array(fringe_table['expnum'], dtype=int)
        fringe_table['ccdnum'] = np.array(fringe_table['ccdnum'], dtype=int)
        # fringe_table["frgscale"] = fringe_table["frgscale"].astype(float)
        # fringe_table['frgscale_median'] = np.nanmedian(fringe_table['frgscale'])
        fringe_table.write(frgscale_output_path, format='ascii.commented_header', overwrite=True)

    print('expnum {} done!'.format(expnum))

    return None


print('=============================================================================')

expnum_list = np.array(exp['expnum'])
np.random.seed(5123490)
expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)

n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(fringe_template_fitting, expnum_list, chunksize=1)
print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))

