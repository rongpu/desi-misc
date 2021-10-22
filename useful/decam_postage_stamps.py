from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

from scipy.ndimage.filters import gaussian_filter

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

ccd_ra = [-0.31244368,-0.00214103, 0.30855858,-0.46789986,-0.1573787 , 0.15336207,
  0.4637642 ,-0.62325889,-0.312972  ,-0.00212455, 0.30866507, 0.61908193,
 -0.77859061,-0.46870955,-0.15780883, 0.15334942, 0.46418217, 0.77441054,
 -0.77876058,-0.46892617,-0.15799484, 0.15333136, 0.46448109, 0.77444204,
 -0.93389515,-0.624237  ,-0.31362077,-0.00213867, 0.30892024, 0.61974856,
  0.92929411,-0.93410772,-0.62439031,-0.31379523,-0.00251046, 0.30860373,
  0.61929563, 0.92907893,-0.77928668,-0.46927775,-0.15819325, 0.15315534,
  0.464108  , 0.77408146,-0.7791703 ,-0.46938561,-0.15825837, 0.15269545,
  0.46382537, 0.77383443,-0.6239286 ,-0.31363566,-0.00262614, 0.30814956,
  0.61848423,-0.46862823,-0.15833137, 0.15254403, 0.46295505,-0.31333245,
  0.30765903]
ccd_dec = [ 0.90299039, 0.90274404, 0.90285652, 0.73894001, 0.73933177, 0.73919444,
  0.73865878, 0.5745655 , 0.57508801, 0.57510357, 0.57486577, 0.57414278,
  0.41001556, 0.41059824, 0.41088721, 0.41057117, 0.41032572, 0.40963196,
  0.24595122, 0.24597951, 0.24624207, 0.24619019, 0.24582139, 0.24534302,
  0.08128957, 0.08150002, 0.08130657, 0.08138846, 0.0810964 , 0.08093379,
  0.08089282,-0.08302691,-0.08319348,-0.08340522,-0.08351659,-0.08366242,
 -0.08355805,-0.08365399,-0.24756494,-0.2479717 ,-0.24812127,-0.24835309,
 -0.2482645 ,-0.2480924 ,-0.41173856,-0.41236738,-0.41281328,-0.41296242,
 -0.41270174,-0.41225407,-0.57638265,-0.57687683,-0.57711492,-0.57725814,
 -0.57674114,-0.74071528,-0.74115162,-0.74130891,-0.74095896,-0.9049206 ,
 -0.90515532]

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

image_dir = '/global/project/projectdirs/cosmo/staging'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
# surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-decam-dr9.fits.gz'
# surveyccd_path_dr8 = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9-garage/reorg/decam/survey-ccds-decam-dr8-newlocs2.fits.gz'
surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr9-trim.fits'
surveyccd_path_dr8 = '/global/cfs/cdirs/desi/users/rongpu/useful/survey-ccds-decam-dr8-trim.fits'

image_vrange = {'u':5, 'g':5, 'r':6, 'i':10, 'z':30, 'Y':30}

################################################################################

def decam_plot(exposure, plot_path, figsize=(13, 12), vrange=None, dr8=False, binsize=20, median=True,
    blob_mask=False, ood_mask=False, gaussian_sigma=None, show=False):
    '''
    Create high-resolution DECam images.
    Example:
    decam_plot(781475, 'tmp_mask_median.jpeg', binsize=20, blob_mask=True, ood_mask=True, median=True) 
    '''
    if (type(exposure)==str) or (type(exposure)==np.str_):
        image_path = os.path.join(image_dir, exposure.strip())
    elif isinstance(exposure, int) or isinstance(exposure, np.integer):
        if not dr8:
            ccd = Table.read(surveyccd_path)
        else:
            ccd = Table.read(surveyccd_path_dr8)
        ccd_index = np.where(ccd['expnum']==exposure)[0][0]
        # band = ccd['image_filename'][ccd_index][ccd['image_filename'][ccd_index].find('_ooi_')+5]
        image_path = os.path.join(image_dir, ccd['image_filename'][ccd_index].strip())

    print(image_path)

    if not os.path.isfile(image_path):
        print('File does not exist:', image_path)
        return None

    if vrange is None:
        band = image_path[image_path.find('_ooi_')+5]
        vrange = image_vrange[band]

    if blob_mask:
        str_loc = str.find(ccd['image_filename'][ccd_index].strip(), '.fits')
        img_filename_base = ccd['image_filename'][ccd_index].strip()[:str_loc]
        blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
        if (not os.path.isfile(blob_path)) or (os.stat(blob_path).st_size==0):
            print(blob_path+' is empty or does not exist!')
            return None
        blob_data = np.load(blob_path)

    plt.figure(figsize=figsize)

    for ii, ccdnum in enumerate(ccdnum_list):

        ccdname = ccdnamenumdict_inv[ccdnum]

        try:
            img = fitsio.read(image_path, ext=ccdname)
        except OSError:
            if ccdname!='S30': # mute S30
                print('{} does not exist in image ({})!'.format(ccdname, exposure))
            continue

        if ood_mask:
            ood_path = image_path.replace('_ooi_', '_ood_')
            ood = fitsio.read(ood_path, ext=ccdname)

        if blob_mask:
            try:
                with fitsio.FITS(image_path) as f:
                    hdu_index = f.movnam_ext(ccdname)
                blob = blob_data['hdu'+str(hdu_index).zfill(2)]
            except (OSError, KeyError):
                print(blob_path+' hdu'+str(hdu_index)+' does not exist!')
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

        # to ignore NAN values, use np.nanmean or np.nanmedian
        if not median:
            img = np.nanmean(np.nanmean(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)
        else:
            img = np.nanmedian(np.nanmedian(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)

        img[~np.isfinite(img)] = 0

        ################################################

        if gaussian_sigma is not None:
            img = gaussian_filter(img, gaussian_sigma, mode='reflect', truncate=3)

        ysize, xsize = img.shape
        ra, dec = ccd_ra[ii], ccd_dec[ii]
        
        fig = plt.imshow(img.T, cmap='seismic', vmin=-vrange, vmax=vrange, 
                   extent=(ra-ysize*pix_size*binsize/2, ra+ysize*pix_size*binsize/2, dec-xsize*pix_size*binsize/2, dec+xsize*pix_size*binsize/2))

    plt.axis([1.07, -1.07, -1.0, 1.0])
    plt.axis('off')
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    # plt.colorbar(fraction=0.04, pad=0.04)
    plt.tight_layout()
    plt.savefig(plot_path)
    if show:
        plt.show()
    else:
        plt.close()


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


def decam_postage_stamp(exposure, binsize=120, plot_path=None, save_path=None, vrange=None, dr8=False, median=True,
    blob_mask=True, ood_mask=True, show=False):
    '''
    Create low-resolution postage stamps.
    Examples:
    decam_postage_stamp(781475, plot_path='tmp_stamp.png', save_path='stamp.npz')
    decam_postage_stamp('/global/cfs/cdirs/cosmo/staging/decam/CP/V4.8.2a/CP20181006/c4d_181007_074137_ooi_g_ls9.fits.fz', plot_path='tmp_stamp.png', save_path='stamp.npz')
    '''
    if type(exposure)==str:
        image_path = exposure
    elif isinstance(exposure, int) or isinstance(exposure, np.integer):
        if not dr8:
            ccd = Table.read(surveyccd_path)
        else:
            ccd = Table.read(surveyccd_path_dr8)
        ccd_index = np.where(ccd['expnum']==exposure)[0][0]
        image_path = os.path.join(image_dir, ccd['image_filename'][ccd_index].strip())
    else:
        raise ValueError('exposure can either be string or integer!')

    print(image_path)
    band = image_path[image_path.find('_ooi_')+5]

    if save_path is not None and os.path.isfile(save_path):

        fullimg = np.load(save_path)['data']

    else:

        if not os.path.isfile(image_path):
            print('File does not exist:', image_path)
            return None

        if blob_mask:
            str_loc = str.find(image_path, '.fits')
            img_filename_base = image_path[len(image_dir)+1:str_loc]
            blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
            try:
                blob_data = np.load(blob_path)
            except FileNotFoundError:
                print(blob_path+' does not exist!')
                return None

        x_pix_small, y_pix_small = np.rint(x_pix/binsize).astype(int), np.rint(y_pix/binsize).astype(int)
        x_size_small, y_size_small = int(np.ceil(x_size/binsize)), int(np.ceil(y_size/binsize))
        fullimg = np.zeros((y_size_small, x_size_small))
        fullimg[:] = np.nan

        for ii, ccdnum in enumerate(ccdnum_list):

            ccdname = ccdnamenumdict_inv[ccdnum]

            try:
                img = fitsio.read(image_path, ext=ccdname)
            except OSError:
                print(ccdname+' does not exist in image!')
                continue

            if ood_mask:
                ood_path = image_path.replace('_ooi_', '_ood_')
                ood = fitsio.read(ood_path, ext=ccdname)

            if blob_mask:
                try:
                    with fitsio.FITS(image_path) as f:
                        hdu_index = f.movnam_ext(ccdname)
                    blob = blob_data['hdu'+str(hdu_index).zfill(2)]
                except (OSError, KeyError):
                    print(blob_path+' hdu'+str(hdu_index)+' does not exist!')
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
            mask = nanmask<0.01 # require at least 1% of the pixels to be unmasked
            img[mask] = np.nan

            ################################################
            img = img.T
            img_y_size, img_x_size = img.shape
            fullimg[y_pix_small[ii]:y_pix_small[ii]+img_y_size, x_pix_small[ii]:x_pix_small[ii]+img_x_size] = img

        fullimg = np.flip(fullimg, 1)

        if save_path is not None:
            np.savez_compressed(save_path, data=fullimg)

    fullimg[~np.isfinite(fullimg)] = 0

    if plot_path is not None:

        if vrange is None:
            vrange = image_vrange[band]

        ax = create_image(fullimg, cmap='seismic', vmin=-vrange, vmax=vrange)
        plt.savefig(plot_path)
        if show:
            plt.show()
        else:
            plt.close()
