# Finish the remaining exposures in the south now that the southern blob masks are available

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits
import healpy as hp
from astropy import wcs

from scipy.ndimage.filters import gaussian_filter
from pathlib import Path
from scipy import stats
from multiprocessing import Pool
import argparse

################################################################################

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
                  'N29': 60, 'N30': 61, 'N31': 62}
ccdnamenumdict_inv = {aa: bb for bb, aa in ccdnamenumdict.items()}

ccdnum_list = [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
               18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
               35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
               52, 53, 54, 55, 56, 57, 58, 59, 60, 62]

img_shape = (4094, 2046)

################################################################################

n_processes = 64

overwrite = False

################################################################################
debug = False
################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('run')
args = parser.parse_args()
run = int(args.run)

blob_dir_dr9 = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
blob_dir_dr10 = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_ccd_blob_mask'

image_dir = '/global/project/projectdirs/cosmo/staging'
##################################################################################################
image_dir_wihtout_cp_fringe = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/undo_cp_fringe_corr'
##################################################################################################

# surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-decam-dr10-v2.fits'
surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey_ccds_in_runs/survey-ccds-decam-dr10-v2-run_{}.fits'.format(run)

template_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_templates'
skyscale_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_scales_tmp/'

skyrun = Table(fitsio.read('/global/cfs/cdirs/desi/users/schlafly/decals/skyrunsdr10-v3.fits'))
print('skyrun', len(skyrun))

mask = skyrun['filter']=='z'
skyrun = skyrun[mask]
print('skyrun', len(skyrun))

exp = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey-ccds-decam-dr10-v2-unique_exps_z_fringe.fits'))
print(len(exp))
skyrun = join(skyrun, exp[['expnum', 'old_fringe', 'new_fringe']])
print('skyrun', len(skyrun))

# # Remove images with CP-only fringe correction
# old_fringe = skyrun['old_fringe']
# new_fringe = skyrun['new_fringe']
# mask = old_fringe & (~new_fringe)
# skyrun = skyrun[~mask]
# print('skyrun', len(skyrun))

mask = skyrun['run']==run
if np.sum(mask)==0:
    sys.exit()
skyrun = skyrun[mask]
print('skyrun', len(skyrun))

expnum_list = skyrun['expnum'].copy()

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)

if debug:
    expnum_list = np.random.choice(expnum_list, size=n_processes, replace=False)

print('Number of exposures in this run:', len(expnum_list))

# ccd_columns = ['image_hdu', 'expnum', 'ccdname', 'ccdskycounts']
ccd = Table(fitsio.read(surveyccd_path))

binsize = 4
pix_size = 0.262/3600*binsize

image_vrange = {'g':5, 'r':6, 'z':30}


# Get run info
mask = skyrun['run']==run
band = skyrun['filter'][mask][0]
run = skyrun['run'][mask][0]
sky_path = os.path.join(template_dir, 'sky_template_{}_{}.fits'.format(band, run))
vrange = image_vrange[band]

sky_dict = {}


print('loading sky template')
for ii, ccdnum in enumerate(ccdnum_list):

    ccdname = ccdnamenumdict_inv[ccdnum]

    if ccdname=='S7' or ccdname=='S30':
        continue

    try:
        sky_dict[ccdname] = fits.getdata(sky_path, extname=ccdname)
    except:
        print(ccdname+' does not exist in template!')
        continue


print('Start!')

def template_fitting(expnum, diagnostic_touch=True):

    start = time.time()

    # print('Expnum', expnum)

    # Get run info
    skyrun_index = np.where(skyrun['expnum']==expnum)[0][0]

    # mask = skyrun['run']==run
    # skyrun_idx = np.where(mask)[0]
    # print('\nrun {}, {} exposures'.format(run, len(skyrun_idx)))

    image_filename = skyrun['image_filename'][skyrun_index].strip()
    image_path = os.path.join(image_dir, image_filename)
    # ood_path = image_path.replace('_ooi_', '_ood_')

    str_loc = str.find(image_filename, '.fits')
    img_filename_base = image_filename[:str_loc]
    skyscale_path = os.path.join(skyscale_dir, image_filename.replace('.fits.fz', '-skyscale.txt'))

    blob_path_dr9 = os.path.join(blob_dir_dr9, 'blob_mask', img_filename_base+'-blobmask.npz')
    blob_path_dr10 = os.path.join(blob_dir_dr10, 'blob_mask', img_filename_base+'-blobmask.npz')

    if os.path.isfile(blob_path_dr9) or (not os.path.isfile(blob_path_dr10)):
        return None

    if (overwrite is False) and os.path.isfile(skyscale_path):
        time_modified = os.path.getmtime(blob_path_dr10)
        if time_modified<1644107000:  # Only rerun fitting if the blob mask is updated
            return None

    ##########################################################################################################
    header_keys = fitsio.read_header(image_path, ext=1).keys()
    if ('FRGSCALE' in header_keys) and ('FRGSCNEW' not in header_keys):
        image_path = os.path.join(image_dir_wihtout_cp_fringe, skyrun['image_filename'][skyrun_index].strip())
    ##########################################################################################################

    if not os.path.exists(os.path.dirname(skyscale_path)):
        try:  # in case another process is also creating the directory
            os.makedirs(os.path.dirname(skyscale_path))
        except:
            pass

    result = Table(names=('image_hdu', 'ccdname', 'ccdskyscale', 'medianskyscale'), dtype=('i4', 'S3', 'f4', 'f4'))
    result['ccdskyscale'].format = '.5f'
    result['medianskyscale'].format = '.5f'

    if os.path.isfile(blob_path_dr9) and (os.stat(blob_path_dr9).st_size!=0):
        blob_path = blob_path_dr9
    elif os.path.isfile(blob_path_dr10) and (os.stat(blob_path_dr10).st_size!=0):
        blob_path = blob_path_dr10
    else:
        print('no blobmask!', skyscale_path)
        if diagnostic_touch:
            Path('/global/u2/r/rongpu/temp/sky_scale_being_written/expnum_{}'.format(expnum)).touch()
        result.write(skyscale_path, format='ascii.commented_header', overwrite=True)
        if diagnostic_touch:
            os.remove('/global/u2/r/rongpu/temp/sky_scale_being_written/expnum_{}'.format(expnum))
        return None

    # try:
    #     blob_data = np.load(blob_path)
    # except:
    #     print(blob_path+' does not exist!')
    #     result.write(skyscale_path, format='ascii.commented_header', overwrite=True)
    #     continue

    blob_data = np.load(blob_path)

    for ii, ccdnum in enumerate(ccdnum_list):

        ccdname = ccdnamenumdict_inv[ccdnum]
        # print(ii, ccdname)

        if ccdname=='S7' or ccdname=='S30':
            continue

        if ccdname not in sky_dict.keys():
            continue

        try:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore')
                img = fits.getdata(image_path, extname=ccdname)
                # ood = fits.getdata(ood_path, extname=ccdname)
        except KeyError:
            print(ccdname+' does not exist in image!')
            continue
        except:
            print('Corrupted?', image_path, expnum, ccdname)
            continue

        sky = sky_dict[ccdname]

        if np.all(sky==0):
            print(ccdname+' template is all zeros!')
            continue

        tmp = np.where((ccd['expnum']==expnum) & (ccd['ccdname']==ccdname))[0]
        if len(tmp)==0:  # no idea why it happens
            print('why why why', expnum, ccdname, image_path)
            continue
        ccd_index = tmp[0]

        # Get HDU index
        # with fitsio.FITS(img_fn) as f:
        #     hdu_index = f.movnam_ext(ccdname)
        hdu_index = ccd['image_hdu'][ccd_index]

        # ######################################################################
        # # Rescale the sky template by ccdskycounts
        # sky *= ccd['ccdskycounts'][ccd_index]
        # ######################################################################

        try:
            blob = blob_data['hdu'+str(hdu_index).zfill(2)]
        except:
            # print(blob_path+' hdu'+str(hdu_index)+' does not exist!')
            continue

        # # Remove median sky
        # sky = np.median(img[blob].flatten())
        # img = img - sky

        # # Apply blob and ood mask
        # img_mask = (blob==True) & (ood==0)
        # Apply blob mask
        img_mask = (blob==True)

        # sky estimation
        median_sky = np.median(img[img_mask].flatten())
        img = img - median_sky

        img1 = img.copy()
        img1[~img_mask] = np.nan
        img_mask = np.isfinite(img1)

        if np.sum(~img_mask)/np.prod(img_mask.shape)>0.9:
            print('{} {:.3f}% pixels are masked; skip'.format(ccdname, np.sum(~img_mask)/np.prod(img_mask.shape)*100))
            continue
        # elif np.sum(~img_mask)/np.prod(img_mask.shape)>0.3:
        #     print('{} {:.3f}% pixels are masked'.format(ccdname, np.sum(~img_mask)/np.prod(img_mask.shape)*100))

        # 3-sigma clipping
        img1_nmad = nmad(img1[img_mask]) # sky level
        # mask = (img1<-3*img1_nmad) | (img1>3*img1_nmad)
        # img1[mask] = 0
        ##############################################################
        mask = (img1<-3*img1_nmad) | (img1>3*img1_nmad)
        # if np.sum(mask)/np.sum(img_mask)>0.03:
        #     print('{} {:.3f}% pixels are clipped'.format(ccdname, np.sum(mask)/np.sum(img_mask)*100))
        ##############################################################
        mask = (img1<-3*img1_nmad)
        img1[mask] = -3*img1_nmad
        mask = (img1>3*img1_nmad)
        img1[mask] = 3*img1_nmad

        img1_flat = img1[img_mask].flatten()
        sky_flat = sky[img_mask].flatten()
        slope, intercept, r_value, p_value, std_err = stats.linregress(sky_flat, img1_flat)
        result.add_row((hdu_index, ccdname, slope, 0))

    n_ccd = len(result)
    if n_ccd>0:
        medianskyscale = np.median(result['ccdskyscale'])
        result['medianskyscale'] = medianskyscale
    else:
        print('No CCD available!')

    if diagnostic_touch:
        Path('/global/u2/r/rongpu/temp/sky_scale_being_written/expnum_{}'.format(expnum)).touch()
    result.write(skyscale_path, format='ascii.commented_header', overwrite=True)
    if diagnostic_touch:
        os.remove('/global/u2/r/rongpu/temp/sky_scale_being_written/expnum_{}'.format(expnum))

    end = time.time()
    print('Expnum {} took {:.1f} seconds'.format(expnum, end - start))


def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(template_fitting, expnum_list, chunksize=1)
        # res = pool.map(template_fitting, expnum_list)

    print('run {} done!!!!!!!!!!!!!!!!!!!!!'.format(run))

if __name__=="__main__":
    main()
