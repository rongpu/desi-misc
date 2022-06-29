# Compute the final fringe templates

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
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

# ccdnum = 5

# # skip S7
# if ccdnum==31:
#     print('skipping S7')
#     sys.exit()

# # skip  S30
# if ccdnum==2:
#     print('skipping S30')
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
# surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/survey-ccds-z-fringe.fits'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_ccd_blob_mask'

plot_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/plots'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/fringe_templates_raw'

output_path = os.path.join(output_dir, 'DECam_z_frg_{}_CCD{}_raw.fits'.format(ccdname, str(ccdnum).zfill(2)))

if os.path.isfile(output_path):
    print('Error: {} already exist!!!'.format(output_path))
    sys.exit()

# Load DR9 fringe template
dr9_fringe_template_path = os.path.join(fringe_dr9_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
if ccdname!='S7':
    fringe = fitsio.read(dr9_fringe_template_path)
    # remove the edge pixels
    fringe = fringe[1:4095, 1:2047]

############################ Load CCD list and weights ############################
ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/survey-ccds-z-fringe.fits'))
print(len(ccd), len(np.unique(ccd['expnum'])))

frgscales = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/frgscale_init/survey-ccds-dr10-deep-fields-v1-z-fringe-stacking-frgscales.fits'))
print(len(frgscales), len(np.unique(frgscales['expnum'])))
mask = np.isfinite(frgscales['slope'])
frgscales = frgscales[mask]
print(len(frgscales), len(np.unique(frgscales['expnum'])))

frgscales['slope_median'] = 0.
for expnum in np.unique(frgscales['expnum']):
    mask = frgscales['expnum']==expnum
    frgscales['slope_median'][mask] = np.nanmedian(frgscales['slope'][mask])
    if not np.all(np.isfinite(frgscales['slope'][mask])):
        print(expnum, np.sum(np.isfinite(frgscales['slope'][mask])), np.sum(mask))

tmp = Table()
tmp['expnum'], tmp['count'] = np.unique(frgscales['expnum'], return_counts=True)
mask = tmp['count']>=55
print(np.sum(~mask))
tmp = tmp[mask]
mask = np.in1d(frgscales['expnum'], tmp['expnum'])
frgscales = frgscales[mask]
print(len(frgscales), len(np.unique(frgscales['expnum'])))
print()

_, idx = np.unique(frgscales['expnum'], return_index=True)
frgscales = frgscales[idx].copy()
print(len(ccd), len(np.unique(ccd['expnum'])))
ccd = join(ccd, frgscales[['expnum', 'slope_median']], keys='expnum', join_type='inner')
print(len(ccd), len(np.unique(ccd['expnum'])))

####################################################################################

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx]
str_loc = np.char.find(np.array(exp['image_filename'], dtype='str'), '/CP20')
exp['obs_date'] = np.array([exp['image_filename'][i][str_loc[i]+3:str_loc[i]+11] for i in range(len(exp))])
print('Total number of nights:', len(np.unique(exp['obs_date'])))

mask = ccd['ccdname']==ccdname
ccd = ccd[mask]
print('{}   Number of images: {}'.format(ccdname, len(ccd)))

# np.random.seed(23152)
# idx = np.random.choice(len(ccd), 2000, replace=False)
# ccd = ccd[idx]
# ccd.sort('exptime')  # Select the longer exposures
# ccd = ccd[-2000:]
# print('{}   Number of images: {}'.format(ccdname, len(ccd)))


def process_ccd(ccd_index):

    if ccd_index%10==0:
        print('{} / {}'.format(ccd_index, len(ccd)))

    hdu_index = ccd['image_hdu'][ccd_index]

    # Load CCD image
    img_fn = os.path.join(image_dir, ccd['image_filename'][ccd_index]).strip()
    ood_fn = img_fn.replace('_ooi_', '_ood_')

    header = fitsio.read_header(img_fn, ext=hdu_index)
    if 'FRGSCNEW' in header.keys():
        frgscale_dr9 = header['FRGSCNEW']
    else:
        # print('FRGSCNEW does not exist!!!')
        frgscale_dr9 = 0
    # w = wcs.WCS(hdulist[ccd['image_hdu'][ccd_index]].header)

    if ('FRGSCALE' in header.keys()) and ('FRGSCNEW' not in header.keys()):
        raise ValueError

    img = fitsio.read(img_fn, ext=hdu_index)
    ood = fitsio.read(ood_fn, ext=hdu_index)

    # Back out the exisiting fringe correction
    if ccdname!='S7':
        img += fringe*frgscale_dr9

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

    # Apply mask
    img[~mask_good] = np.nan

    # Only keep the good half of the S7
    if ccdname=='S7':
        half = img_shape[1] // 2
        img = img[:, :half]

    # Remove median sky
    sky = np.nanmedian(img.flatten())
    img = img - sky

    # Remove the constant offset between the amps
    if ccdname!='S7':
        img1 = img[:, img_shape[1]//2-10:img_shape[1]//2]
        if np.sum(np.isfinite(img1))>100:
            img[:, :img_shape[1]//2] -= np.nanmedian(img1)
        img1 = img[:, img_shape[1]//2:img_shape[1]//2+10]
        if np.sum(np.isfinite(img1))>100:
            img[:, img_shape[1]//2:] -= np.nanmedian(img1)

    # Normalize by initial frgscale
    frgscale = ccd['slope_median'][ccd_index]
    img = img/frgscale

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

    # trim edges
    if ccdname!='S7':
        trim_size_x, trim_size_y = 15, 15   # trimmed imagee size need to be multiples of binsize
        img1 = img1[trim_size_y:(img1.shape[0]-trim_size_y), trim_size_x:(img1.shape[1]-trim_size_x)]
    else:
        trim_size_x, trim_size_y = 0, 1
        img1 = img1[trim_size_y:(img1.shape[0]-trim_size_y), trim_size_x:(img1.shape[1]-trim_size_x)]
    # downsize the image to speed up gaussian filter
    if ccdname!='S7':
        binsize = 32
    else:
        binsize = 31
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        img1 = np.nanmedian(np.nanmedian(img1.reshape((img1.shape[0]//binsize, binsize, img1.shape[1]//binsize,-1)), axis=3), axis=1)
    x_small_grid = trim_size_x + binsize/2+binsize*np.arange(img1.shape[1])
    y_small_grid = trim_size_y + binsize/2+binsize*np.arange(img1.shape[0])
    img1[np.isnan(img1)] = 0

    # Gaussian filtering
    img1_smooth = gaussian_filter(img1, 3, mode='reflect')

    # Convert the downsized smooth image to full size
    interp_func = interp2d(x_small_grid, y_small_grid, img1_smooth, kind='linear')
    x_grid, y_grid = np.arange(img.shape[1]), np.arange(img.shape[0])
    smooth_sky = interp_func(x_grid, y_grid).reshape(img.shape)

    # Subtract off the smooth component
    img -= smooth_sky

    # Add back the other half
    if ccdname=='S7':
        tmp = np.zeros(img_shape)
        half = img_shape[1] // 2
        tmp[:, :half] = img
        tmp[:, half:] = np.nan
        img = tmp

    gc.collect()

    return img, frgscale


print('=============================================================================')

n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(process_ccd, np.arange(len(ccd)), chunksize=1)
print('res', len(res))

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)
print('res', len(res))

print('Almost done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))

# img_median_final = np.nanmedian(img_list_all, axis=0)

# Reduce memory usage by chopping up image
def get_median(ii):
    # if ii%(img_shape[0]//50)==0:
    #     print('{}/{}'.format(ii, img_shape[0]))
    tmp = []
    for jj in range(len(res)):
        n_repeat = np.maximum(1, int(res[jj][1]//500))
        for kk in range(n_repeat):
            tmp.append(res[jj][0][ii])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        tmp_median = np.nanmedian(tmp, axis=0)
    gc.collect()
    return tmp_median


n_processes = 64
with Pool(processes=n_processes) as pool:
    res = pool.map(get_median, np.arange(img_shape[0]), chunksize=1)
img_median_final = np.vstack(res)

# Add back the edge pixels
img_median_final = np.pad(img_median_final, ((1, 1), (1, 1)), mode='constant', constant_values=np.nan)

fitsio.write(output_path, img_median_final, clobber=True)


sys.path.append(os.path.expanduser('~/git/Python/useful/'))
from decam_postage_stamps import create_image

vrange = 3e-3
create_image((img_median_final).T, cmap='seismic', vmin=-vrange, vmax=vrange)
plt.savefig(os.path.join(plot_dir, 'DECam_z_frg_{}_CCD{}_raw.png'.format(ccdname, str(ccdnum).zfill(2))))
plt.close()

vrange = 1e-3
if ccdname!='S7':
    create_image((fringe).T, cmap='seismic', vmin=-vrange, vmax=vrange)
    plt.savefig(os.path.join(plot_dir, 'DECam_z_frg_{}_CCD{}_dr9.png'.format(ccdname, str(ccdnum).zfill(2))))
    plt.close()

print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))
print('@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@')
