from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

from pathlib import Path
from multiprocessing import Pool
import argparse

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'} 
plt.rcParams.update(params)

# nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

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

x_size, y_size = 29590, 26787 # Pixel size if all the CCDs are stitched into one image
x_pix = np.array([ 8498. , 12747. , 16996. ,  6373.5, 10622.5, 14871.5, 19120.5,
                   4249. ,  8498. , 12747. , 16996. , 21245. ,  2124.5,  6373.5,
                  10622.5, 14871.5, 19120.5, 23369.5,  2124.5,  6373.5, 10622.5,
                  14871.5, 19120.5, 23369.5,    -0. ,  4249. ,  8498. , 12747. ,
                  16996. , 21245. , 25494. ,    -0. ,  4249. ,  8498. , 12747. ,
                  16996. , 21245. , 25494. ,  2124.5,  6373.5, 10622.5, 14871.5,
                  19120.5, 23369.5,  2124.5,  6373.5, 10622.5, 14871.5, 19120.5,
                  23369.5,  4249. ,  8498. , 12747. , 16996. , 21245. ,  6373.5,
                  10622.5, 14871.5, 19120.5,  8498. , 16996. ])
y_pix = np.array([    0,     0,     0,  2249,  2249,  2249,  2249,  4498,  4498,
                   4498,  4498,  4498,  6747,  6747,  6747,  6747,  6747,  6747,
                   8996,  8996,  8996,  8996,  8996,  8996, 11245, 11245, 11245,
                  11245, 11245, 11245, 11245, 13494, 13494, 13494, 13494, 13494,
                  13494, 13494, 15743, 15743, 15743, 15743, 15743, 15743, 17992,
                  17992, 17992, 17992, 17992, 17992, 20241, 20241, 20241, 20241,
                  20241, 22490, 22490, 22490, 22490, 24739, 24739])

img_shape = (4094, 2046)
pix_size = 0.262/3600

############################################################

n_processes = 32

image_vrange = {'g':0.04, 'r':0.05, 'z':0.3} # nmgy / sq. arcsec.
image_vrange_raw = {'g':3, 'r':4, 'z':20} # raw image values

image_dir = '/global/project/projectdirs/cosmo/staging'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
# surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'
# surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-trim.fits'
surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-trim-less.fits'
# surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9-garage/reorg/decam/survey-ccds-decam-dr8-newlocs2.fits.gz'
# surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr8-trim.fits'
sky_scales_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/calib/sky_pattern/sky-scales.fits'
sky_template_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/calib/sky_pattern/sky_templates'

save_dir = '/global/cscratch1/sd/rongpu/dr10dev/decam_postage/postage_save'
plot_dir = '/global/cscratch1/sd/rongpu/dr10dev/decam_postage/postage_image'
touch_dir = '/global/cscratch1/sd/rongpu/dr10dev/decam_postage/postage_being_written'


parser = argparse.ArgumentParser()
parser.add_argument('args')
args = parser.parse_args()
n_task, task_id = [int(tmp) for tmp in args.args.split()]

ccd = Table.read(surveyccd_path)
skyscales = Table.read(sky_scales_path)

expnum_list = np.sort(np.unique(ccd['expnum']))
np.random.seed(261)
expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)

# ############################################################
# mask_ccd = ccd['ccd_cuts_ok']==True
# expnum_list = np.unique(ccd['expnum'][mask_ccd])
# np.random.seed(261)
# expnum_list = np.random.choice(expnum_list, size=1024, replace=False)
# ############################################################

# split among the Cori nodes
expnum_list_split = np.array_split(expnum_list, n_task)
expnum_list = expnum_list_split[task_id]
print('Number of expousres (task_id {}): {}'.format(task_id, len(expnum_list)))


def create_image(data, cmap='gray', dpi=80, vmin=None, vmax=None, origin=None, norm=None):
    '''
    Create an image with exactly the same pixel dimension as the data.
    Example:
        x, y = np.arange(0, 10), np.arange(0, 10)
        xx, yy = np.meshgrid(x, y)
        img = np.array((xx + yy)%2==0, dtype=int)
        ax = create_image(img)
        plt.savefig('img.png')
        plt.close()
    '''
    xpixels, ypixels = data.shape[0], data.shape[1]
    figsize = ypixels / dpi, xpixels / dpi
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(data, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, origin=origin, norm=norm)
    plt.axis('off')
    return ax


def decam_postage_stamp(expnum, binsize=120, save_path=None, dr8=False, median=True,
    blob_mask=True, ood_mask=True, show=False, sky_correction=True, diagnostic_touch=True):
    '''
    Examples:
    decam_postage_stamp(781475, save_path='stamp.npz')
    '''
    ccd_index = np.where(ccd['expnum']==expnum)[0][0]
    image_filename = ccd['image_filename'][ccd_index].strip()

    image_path = os.path.join(image_dir, image_filename)
    # print(image_path)
    band = image_path[image_path.find('_ooi_')+5]

    save_path = os.path.join(save_dir, image_filename.replace('.fits.fz', '.npz'))
    plot_path = os.path.join(plot_dir, image_filename.replace('.fits.fz', '.png'))
    if not os.path.exists(os.path.dirname(save_path)):
        try:
            os.makedirs(os.path.dirname(save_path))
        except:
            pass
    if not os.path.exists(os.path.dirname(plot_path)):
        try:
            os.makedirs(os.path.dirname(plot_path))
        except:
            pass

    if os.path.isfile(save_path):

        fullimg = np.load(save_path)['data']

    else:

        if not os.path.isfile(image_path):
            print('File does not exist:', image_path)
            return None

        if blob_mask:
            str_loc = str.find(image_path, '.fits')
            img_filename_base = image_path[len(image_dir)+1:str_loc]
            blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
            if os.path.isfile(blob_path) and (os.stat(blob_path).st_size>0):
                try:
                    blob_data = np.load(blob_path)
                except FileNotFoundError:
                    print(blob_path+' does not exist ({})!'.format(expnum))
                    return None
            else:
                blob_mask = False

        x_pix_small, y_pix_small = np.rint(x_pix/binsize).astype(int), np.rint(y_pix/binsize).astype(int)
        x_size_small, y_size_small = int(np.ceil(x_size/binsize)), int(np.ceil(y_size/binsize))
        fullimg = np.zeros((y_size_small, x_size_small))
        fullimg[:] = np.nan

        for ii, ccdnum in enumerate(ccdnum_list):

            # #################
            # print(ii, ccdnum)
            # #################

            ccdname = ccdnamenumdict_inv[ccdnum]

            try:
                img = fitsio.read(image_path, ext=ccdname)
            except OSError:
                if ccdname!='S30':  # mute S30
                    print('{} does not exist in image ({})!'.format(ccdname, expnum))
                continue

            if sky_correction:
                sky_index = np.where((skyscales['expnum']==expnum) & (skyscales['ccdname']==ccdname))[0][0]
                run = skyscales['run'][sky_index]
                skyscale = skyscales['skyscale'][sky_index]
                if (run!=-1) and (skyscale!=0):
                    sky_template_path = os.path.join(sky_template_dir, 'sky_template_{}_{}.fits.fz'.format(band, run))
                    try:
                        sky_template = fitsio.read(sky_template_path, ext=ccdname)
                    except OSError:
                        raise OSError('{} does not exist in sky template {} ({})!'.format(ccdname, run, expnum))
                    img -= skyscale * sky_template

            if ood_mask:
                ood_path = image_path.replace('_ooi_', '_ood_')
                ood = fitsio.read(ood_path, ext=ccdname)

            if blob_mask:
                try:
                    with fitsio.FITS(image_path) as f:
                        hdu_index = f.movnam_ext(ccdname)
                    blob = blob_data['hdu'+str(hdu_index).zfill(2)]
                except (OSError, KeyError):
                    print(blob_path+' hdu'+str(hdu_index)+' does not exist ({})!'.format(expnum))
                    continue

            if ood_mask or blob_mask:
                img_mask = np.ones(img.shape, dtype=bool)
                if ood_mask:
                    img_mask &= (ood==0)
                if blob_mask:
                    img_mask &= (blob==True)
                img[~img_mask] = np.nan

            # Only keep the good half of the S7
            if ccdname=='S7':
                img_original = img.copy()
                half = img_shape[1] // 2
                img = img[:, :half]

            # Remove constant background
            if not blob_mask:
                # naive sky estimation
                mask = (img<np.nanpercentile(img.flatten(), 95))
                median_sky = np.nanmedian(img[mask].flatten())
            else:
                median_sky = np.nanmedian(img)
            img = img - median_sky

            # Add back the other half
            if ccdname=='S7':
                tmp = img_original
                half = img_shape[1] // 2
                tmp[:, :half] = img
                tmp[:, half:] = np.nan
                img = tmp

            ################ downsize image ################

            # trim edges to enable downsizing
            # trimmed image size need to be multiples of binsize
            trim_size_x = img.shape[1] % binsize
            trim_size_y = img.shape[0] % binsize
            img = img[:(img.shape[0]-trim_size_y), :(img.shape[1]-trim_size_x)]
            nanmask = np.isfinite(img)

            # to ignore NAN values, use np.nanmean or np.nanmedian
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                if not median:
                    img = np.nanmean(np.nanmean(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)
                else:
                    img = np.nanmedian(np.nanmedian(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)
                nanmask = np.mean(np.mean(nanmask.reshape((nanmask.shape[0]//binsize, binsize, nanmask.shape[1]//binsize,-1)), axis=3), axis=1)
            mask = nanmask<0.01  # require at least 1% of the pixels to be unmasked
            img[mask] = np.nan

            ################################################

            img = img.T
            img_y_size, img_x_size = img.shape
            fullimg[y_pix_small[ii]:y_pix_small[ii]+img_y_size, x_pix_small[ii]:x_pix_small[ii]+img_x_size] = img

        fullimg = np.flip(fullimg, 1)

        if diagnostic_touch:
            Path(os.path.join(touch_dir, 'expnum_{}'.format(expnum))).touch()
        np.savez_compressed(save_path, data=fullimg)
        if diagnostic_touch:
            os.remove(os.path.join(touch_dir, 'expnum_{}'.format(expnum)))

    fullimg[~np.isfinite(fullimg)] = 0

    # convert to nanomaggie per sq. arcsec
    if ccd['median_ccdzpt'][ccd_index]!=0:
        fullimg = fullimg / 10**((ccd['median_ccdzpt'][ccd_index]-22.5)/2.5) / ccd['exptime'][ccd_index] / (0.262**2)
        vrange = image_vrange[band]
    else:
        vrange = image_vrange_raw[band]

    ax = create_image(fullimg, cmap='seismic', vmin=-vrange, vmax=vrange)
    plt.savefig(plot_path)
    if show:
        plt.show()
    else:
        plt.close()


def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(decam_postage_stamp, expnum_list)

    # make_plots(expnum_list[0])

    print('task_id {} done!!!!!!!!!!!!!!!!!!!!!'.format(task_id))


if __name__=="__main__":
    main()
