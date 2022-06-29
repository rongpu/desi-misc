from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits
from astropy import wcs

from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter

import argparse

from multiprocessing import Pool


time_start = time.time()

parser = argparse.ArgumentParser()
parser.add_argument('ccdnum')
args = parser.parse_args()
ccdnum = int(args.ccdnum)

# ccdnum = 27

# # skip S7
# if hdu_index==31:
#     print('skipping S7')
#     sys.exit()

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

ccdname = ccdnamenumdict_inv[ccdnum]
print('CCD: ', ccdname)

fringe_dr9_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/fringe/DECam_CP-Fringe-Normed'
# fringe_cp_dir = '/global/homes/d/djschleg/cosmo/staging/decam/DECam_CP-Fringe'
image_dir = '/global/cfs/cdirs/cosmo/staging/'
surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/misc/survey-ccds-dr10-deep-fields-v1-z-fringe-stacking.fits'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'

dr9_fringe_scale_path = '/global/cfs/cdirs/desi/users/rongpu/data/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-dr9-frgscale.fits'

plot_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/plots'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data'

output_path = os.path.join(output_dir, 'DECam_z_frg_{}_CCD{}_raw.fits'.format(ccdname, str(ccdnum).zfill(2)))

# if os.path.isfile(output_path):
#     print('Error: {} already exist!!!'.format(output_path))
#     sys.exit()

# Load DR9 fringe template
dr9_fringe_template_path = os.path.join(fringe_dr9_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
fringe = fitsio.read(dr9_fringe_template_path)
# remove the edge pixels
fringe = fringe[1:4095, 1:2047]

dr9_frgscale = Table(fitsio.read(dr9_fringe_scale_path))
print('frgscale', len(dr9_frgscale))
# _, idx = np.unique(dr9_frgscale['expnum'], return_index=True)
# dr9_frgscale = dr9_frgscale[idx]
mask = dr9_frgscale['frgscale_new']!=-99.
dr9_frgscale = dr9_frgscale[mask]
print('frgscale', len(dr9_frgscale))
mask = dr9_frgscale['n_fringe_new']>55
dr9_frgscale = dr9_frgscale[mask]
print('frgscale', len(dr9_frgscale))
dr9_frgscale['ccd_id_str'] = np.char.add(np.array(dr9_frgscale['expnum']).astype(str), dr9_frgscale['ccdname'])

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'ccdname', 'filter', 'expnum', 'exptime']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
print(len(ccd))

ccd['ccd_id_str'] = np.char.add(np.array(ccd['expnum']).astype(str), ccd['ccdname'])
mask = np.in1d(ccd['ccd_id_str'], dr9_frgscale['ccd_id_str'])
ccd = ccd[mask]
print(len(ccd))
ccd = join(ccd, dr9_frgscale[['ccd_id_str', 'frgscale_new', 'frgscale_new_median']], keys='ccd_id_str', join_type='left')
print(len(ccd))

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx]
str_loc = np.char.find(np.array(exp['image_filename'], dtype='str'), '/CP20')
exp['obs_date'] = np.array([exp['image_filename'][i][str_loc[i]+3:str_loc[i]+11] for i in range(len(exp))])
print('Total number of nights:', len(np.unique(exp['obs_date'])))

# img_list_all = []

# # skip S7
# if hdu_index==31:
#     sys.exit()

mask = ccd['ccdname']==ccdname
ccd = ccd[mask]
print('{}   Number of images: {}'.format(ccdname, len(ccd)))

# idx = np.random.choice(len(ccd), 2000, replace=False)
# ccd = ccd[idx]

# for ccd_index in range(len(ccd)):
def process_ccd(ccd_index):

    if ccd_index%10==0:
        print('{} / {}'.format(ccd_index, len(ccd)))

    hdu_index = ccd['image_hdu'][ccd_index]
    frgscale1 = ccd['frgscale_new'][ccd_index]
    frgscale_median = ccd['frgscale_new_median'][ccd_index]

    # Load CCD image
    img_fn = os.path.join(image_dir, ccd['image_filename'][ccd_index]).strip()
    ood_fn = img_fn.replace('_ooi_', '_ood_')

    header = fitsio.read_header(img_fn, ext=hdu_index)
    if 'FRGSCNEW' in header.keys():
        frgscale = header['FRGSCNEW']
        if frgscale1!=frgscale:
            print(frgscale1, frgscale)
            raise ValueError
    else:
        print('FRGSCNEW does not exist!!!')
        frgscale = 0
    # w = wcs.WCS(hdulist[ccd['image_hdu'][ccd_index]].header)

    img = fitsio.read(img_fn, ext=hdu_index)
    ood = fitsio.read(ood_fn, ext=hdu_index)

    # Back out the exisiting fringe correction
    img += fringe*frgscale

    # Load blob mask
    str_loc = str.find(ccd['image_filename'][ccd_index].strip(), '.fits')
    img_filename_base = ccd['image_filename'][ccd_index].strip()[:str_loc]
    blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
    try:
        blob_data = np.load(blob_path)
    except:
        print('Error loading {}'.format(blob_path))
        return None
    try:
        blob = blob_data['hdu'+str(hdu_index).zfill(2)]
    except:
        print('Error loading {} blobmask'.format(ccdname))
        return None

    mask_good = blob & (ood==0)

    # Remove median sky
    sky = np.median(img[mask_good].flatten())
    img = img - sky

    # Normalize by frgscale
    img = img/frgscale_median

    # Apply mask
    img[~mask_good] = np.nan

    if np.all(np.isnan(img)):
        print('{} CCD {} are all NaNs'.format(img_fn, ccdname))

    # 5-sigma clipping
    sky_nmad = nmad(img)
    mask = (img<-5*sky_nmad)
    img[mask] = -5*sky_nmad
    mask = (img>5*sky_nmad)
    img[mask] = 5*sky_nmad

    #################################### Subtract smooth sky ####################################

    img1 = img.copy()

    # 5-sigma clipping
    sky_nmad = nmad(img1)
    mask = (img1<-5*sky_nmad)
    img1[mask] = -5*sky_nmad
    mask = (img1>5*sky_nmad)
    img1[mask] = 5*sky_nmad

    # trim edges
    trim_size = 11   # trimmed imagee size need to be multiples of binsize
    img1 = img1[trim_size:(img1.shape[0]-trim_size), trim_size:(img1.shape[1]-trim_size)]

    # downsize the image to speed up gaussian filter
    binsize = 4
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        img1 = np.nanmean(np.nanmean(img1.reshape((img1.shape[0]//binsize, binsize, img1.shape[1]//binsize,-1)), axis=3), axis=1)
    x_small_grid = trim_size + binsize/2+binsize*np.arange(img1.shape[1])
    y_small_grid = trim_size + binsize/2+binsize*np.arange(img1.shape[0])
    img1[np.isnan(img1)] = 0

    # Gaussian filtering
    img1_smooth = gaussian_filter(img1, 30, mode='reflect')

    # Convert the downsized smooth image to full size
    interp_func = interp2d(x_small_grid, y_small_grid, img1_smooth, kind='linear')
    x_grid, y_grid = np.arange(img.shape[1]), np.arange(img.shape[0])
    smooth_sky = interp_func(x_grid, y_grid).reshape(img.shape)

    # Subtract off the smooth component
    img -= smooth_sky

    gc.collect()

    return img


print('=============================================================================')

n_processes = 128
with Pool(processes=n_processes) as pool:
    img_list_all = pool.map(process_ccd, np.arange(len(ccd)), chunksize=1)
print('img_list_all', len(img_list_all))

# Remove None elements from the list
for index in range(len(img_list_all)-1, -1, -1):
    if img_list_all[index] is None:
        img_list_all.pop(index)
print('img_list_all', len(img_list_all))

print('Almost done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))

vrange = 4e-3

# img_median_final = np.nanmedian(img_list_all, axis=0)

# Reduce memory usage by chopping up image
def get_median(ii):
    # if ii%(img_shape[0]//50)==0:
    #     print('{}/{}'.format(ii, img_shape[0]))
    tmp = []
    for jj in range(len(img_list_all)):
        tmp.append(img_list_all[jj][ii])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tmp_median = np.nanmedian(tmp, axis=0)
    gc.collect()
    return tmp_median


n_processes = 64
with Pool(processes=n_processes) as pool:
    res = pool.map(get_median, np.arange(img_shape[0]), chunksize=1)
img_median_final = np.vstack(res)

fitsio.write(output_path, img_median_final, clobber=True)


sys.path.append(os.path.expanduser('~/git/Python/useful/'))
from decam_postage_stamps import create_image

create_image((img_median_final).T, cmap='seismic', vmin=-vrange, vmax=vrange)
plt.savefig(os.path.join(plot_dir, 'fringe_new_{}.png'.format(ccdnum)))
plt.close()

# Plot 3-pixel gaussian smoothed fringe image
# 5-sigma clipping
sky_nmad = nmad(img_median_final[np.isfinite(img_median_final)])  # sky level
img_median_final1 = img_median_final.copy()
# mask = (img_median_final<-5*sky_nmad) | (img_median_final>5*sky_nmad)
# img_median_final1[mask] = 0
mask = (img_median_final<-5*sky_nmad)
img_median_final1[mask] = -5*sky_nmad
mask = (img_median_final>5*sky_nmad)
img_median_final1[mask] = 5*sky_nmad
img_median_final1[np.isnan(img_median_final1)] = 0
img_median_3pix_gauss = gaussian_filter((img_median_final1), 3, mode='reflect')
create_image((img_median_3pix_gauss).T, cmap='seismic', vmin=-vrange, vmax=vrange)
plt.savefig(os.path.join(plot_dir, 'fringe_new_{}_smooth.png'.format(ccdnum)))
plt.close()

create_image((fringe).T, cmap='seismic', vmin=-vrange, vmax=vrange)
plt.savefig(os.path.join(plot_dir, 'fringe_old_{}.png'.format(ccdnum)))
plt.close()

print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))
print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
