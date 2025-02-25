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

from pathlib import Path
from multiprocessing import Pool

from scipy import stats

################################################################################

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

img_shape = (4094, 2046)

################################################################################

ccdnum_edge_list = [1, 3, 4, 7, 8, 12, 13, 18, 19, 24, 25, 32, 38, 39, 44, 45, 50, 51, 55, 56, 59, 60, 62]

n_processes = 12

smoothing_scale = 60 # in pixels

cp_image_dir = '/global/project/projectdirs/cosmo/staging'
lg_image_dir = '/global/cfs/cdirs/cosmo/staging/decam/CP-LG9'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'

template_dir = '/global/cscratch1/sd/rongpu/dr9dev/pupil_pattern/pupil_templates/'

pix_size = 0.262/3600

# image_path_list = glob.glob(os.path.join(cp_image_dir, '*ooi*.fits.fz'))
ccd = Table.read('/global/cfs/cdirs/desi/users/rongpu/dr9/pupil_pattern/survey-ccds-decam-dr9-lg-pupilrun.fits')
mask = ccd['blacklist']==False
ccd = ccd[mask]

run_list = np.unique(ccd['run'])
print(len(run_list))

def compute_smoothed_pupil(run):
    
    mask = ccd['run']==run
    band = ccd['filter'][mask][0]

    print('band: {}, run: {}'.format(band, run))

    raw_path = os.path.join(template_dir, 'pupil_raw_{}_{}.fits.fz'.format(band, run))
    smoothed_path = os.path.join(template_dir, 'pupil_smooth_{}_{}.fits.fz'.format(band, run))

    if os.path.isfile(smoothed_path):
        print(smoothed_path+' already exists!')
        return None

    hdul_smoothed = fitsio.FITS(smoothed_path, mode='rw', clobber=True)
    hdul_smoothed.write(data=None) # first HDU is empty

    for ii, ccdnum in enumerate(ccdnum_list):

        img_list = []
        ccdname = ccdnamenumdict_inv[ccdnum]

        try:
            img = fitsio.read(raw_path, ext=ccdname)
        except:
            print('{} does not exist!'.format(ccdname))
            continue

        if (ccdname=='S7'):
            # Only keep the good half of the S7
            half = img_shape[1] // 2
            img = img[:, :half]

        img_original = img.copy()
        img = img

        # Fill in NAN values
        mask = ~np.isfinite(img)
        # print('number of NAN pixels:', np.sum(mask))
        img[mask] = 0

        # 5-sigma clipping
        img_median_value = np.median(img)
        sky_nmad = nmad(img[np.isfinite(img)]) # sky level
        mask = (img<img_median_value-5*sky_nmad) | (img>img_median_value+5*sky_nmad)
        print('5-sigma clipped pixels: {} ({:.2f}%)'.format(np.sum(mask), np.sum(mask)/len(mask.flatten())*100))
        # img[mask] = 0
        mask = (img<img_median_value-5*sky_nmad)
        img[mask] = img_median_value-5*sky_nmad
        mask = (img>img_median_value+5*sky_nmad)
        img[mask] = img_median_value+5*sky_nmad

        if (ccdname=='S7'):

            # trim three of edges
            trim_size = 20
            trim_size_bottom = 1
            img = img[trim_size:(img.shape[0]-trim_size), trim_size:(img.shape[1]-trim_size_bottom)]

            # downsize the image to speed up gaussian filter
            binsize = 1 # no downsizing
            # img = np.nanmean(np.nanmean(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)
            x_small_grid = trim_size + binsize/2 + binsize*np.arange(img.shape[1])
            y_small_grid = trim_size + binsize/2 + binsize*np.arange(img.shape[0])

        else:

            ################### Normal exposures ####################

            # trim edges
            trim_size = 10            
            img = img[trim_size:(img.shape[0]-trim_size), trim_size:(img.shape[1]-trim_size)]

            # downsize the image to speed up gaussian filter
            binsize = 1 # no downsizing
            # img = np.nanmean(np.nanmean(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)
            x_small_grid = trim_size + binsize/2 + binsize*np.arange(img.shape[1])
            y_small_grid = trim_size + binsize/2 + binsize*np.arange(img.shape[0])

        # Gaussian filtering
        img1_smooth = gaussian_filter(img, smoothing_scale/binsize, mode='reflect')

        # Convert the downsized smooth image to full size
        interp_func = interp2d(x_small_grid, y_small_grid, img1_smooth, kind='linear')
        x_grid, y_grid = np.arange(img_original.shape[1]), np.arange(img_original.shape[0])
        img_smooth = interp_func(x_grid, y_grid).reshape(img_original.shape)

        if (ccdname=='S7'):
            # Add back the other half
            tmp = np.zeros(img_shape)
            half = img_shape[1] // 2
            tmp[:, :half] = img_smooth
            img_smooth = tmp

        ################################ Save sky template ###############################

        hdul_smoothed.write(data=img_smooth, extname=ccdname, compress='rice')

        # ##################
        # end = time.time()
        # print('Took {:.1f} seconds'.format(end - start))
        # ##################

    hdul_smoothed.close()
    
def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(compute_smoothed_pupil, run_list)
    
    # make_plots('c4d_181102_051559_ooi_z_lg9.fits.fz')

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

