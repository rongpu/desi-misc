from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
import healpy as hp
from astropy import wcs

from scipy.ndimage import gaussian_filter
from pathlib import Path
from multiprocessing import Pool
from itertools import repeat
import argparse

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

n_processes = 128

binsize = 4
pix_size = 0.262/3600*binsize

image_vrange = {'g':2.5, 'r':3, 'z':15}

overwrite = False

image_dir = '/dvs_ro/cfs/cdirs/cosmo/staging/'
##################################################################################################
image_dir_wihtout_cp_fringe = '/dvs_ro/cfs/cdirs/desi/users/rongpu/data/dr11dev/undo_cp_fringe_corr'
##################################################################################################

template_dir = '/dvs_ro/cfs/cdirs/desi/users/rongpu/data/dr11dev/sky_templates'

plot_dir = '/global/cfs/cdirs/cosmo/www/temp/rongpu/dr11/sky_correction_before_and_after'


def make_plots(expnum, apply_subtraction=True):

    # Get run info
    ss_index = np.where(ss['expnum']==expnum)[0][0]
    band = ss['filter'][ss_index]
    run = ss['run'][ss_index]
    # vrange = image_vrange[band]

    sky_path = os.path.join(template_dir, 'sky_template_{}_{}.fits'.format(band, run))
    if not os.path.isfile(sky_path):
        sky_path = os.path.join(template_dir, 'sky_template_{}_{}.fits.fz'.format(band, run))

    if apply_subtraction:
        plot_path = os.path.join(plot_dir, '{}_{}_image_{}_after_subtraction.png'.format(band, run, expnum))
    else:
        plot_path = os.path.join(plot_dir, '{}_{}_image_{}_original.png'.format(band, run, expnum))

    if (overwrite==False) and os.path.isfile(plot_path):
        # print(plot_path, 'already exists!!!')
        return None

    # surveyccd_path = '/dvs_ro/cfs/cdirs/desi/users/rongpu/data/dr11dev/sky_pattern/survey_ccds_in_runs/survey-ccds-decam-dr11-merged-incl-early-run_{}.fits'.format(run)
    surveyccd_path = '/pscratch/sd/r/rongpu/tmp/survey_ccds_in_runs_for_qa/survey-ccds-decam-dr11-merged-incl-early-run_{}.fits'.format(run)
    ccd_columns = ['image_filename', 'expnum', 'ccdname', 'ccdskycounts', 'ccd_cuts', 'plver', 'zpt', 'exptime']
    ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
    ccd['plver'] = np.char.strip(ccd['plver'])

    # print(len(ccd))

    ccd_index = np.where((ccd['expnum']==expnum))[0][0]
    image_filename = ccd['image_filename'][ccd_index].strip()
    zpt, exptime = ccd['zpt'][ccd_index], ccd['exptime'][ccd_index]
    vrange = (0.5 * 0.262**2) * (10**((zpt-22.5)/2.5) * exptime)  # normalize to +- 0.5 nanomaggies per arcsec^2

    image_path = os.path.join(image_dir_wihtout_cp_fringe, image_filename)
    if not os.path.isfile(image_path):
        image_path = os.path.join(image_dir, image_filename)

    if not os.path.exists(os.path.dirname(plot_path)):
        try:
            os.makedirs(os.path.dirname(plot_path))
        except:
            pass

    if apply_subtraction:

        skyscale = ss['skyscale'][ss_index]
        template_fits = fitsio.FITS(sky_path)

        sky_dict = {}
        for ii, ccdnum in enumerate(ccdnum_list):
            ccdname = ccdnamenumdict_inv[ccdnum]
            if ccdname in template_fits:
                sky_dict[ccdname] = fitsio.read(sky_path, ext=ccdname)

    mask = ccd['expnum']==ss['expnum'][ss_index]
    n_good_ccd = np.sum(ccd['ccd_cuts'][mask]==0)

    print(plot_path)
    # Path(plot_path).touch()

    image_fits = fitsio.FITS(image_path)

    plt.figure(figsize=(13.7, 13.075))

    for ii, ccdnum in enumerate(ccdnum_list):

        # print(ii)
        ccdname = ccdnamenumdict_inv[ccdnum]

        if ccdname in image_fits:
            img = fitsio.read(image_path, ext=ccdname)
        else:
            continue

        # naive sky estimation
        mask = (img<np.percentile(img.flatten(), 95))
        median_sky = np.median(img[mask].flatten())
        img = img - median_sky

        # Apply sky pattern correction
        if apply_subtraction and (ccdname in sky_dict.keys()):
            img = img - sky_dict[ccdname] * skyscale

        ################ downsize image ################

        # trim edges to enable downsizing
        # trimmed image size need to be multiples of binsize
        trim_size_x = img.shape[1] % binsize
        trim_size_y = img.shape[0] % binsize
        img = img[:(img.shape[0]-trim_size_y), :(img.shape[1]-trim_size_x)]

        # to ignore NAN values, use np.nanmean
        img = np.nanmean(np.nanmean(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)

        ################################################

        img[~np.isfinite(img)] = 0
        img = gaussian_filter(img, 3, mode='reflect', truncate=3)

        ysize, xsize = img.shape
        ra, dec = ccd_ra[ii], ccd_dec[ii]

        fig = plt.imshow(img.T, cmap='gray', vmin=-vrange, vmax=vrange,
                   extent=(ra-ysize*pix_size/2, ra+ysize*pix_size/2, dec-xsize*pix_size/2, dec+xsize*pix_size/2))

    text = '{} band, expnum = {}\n'.format(band, expnum)
    text += 'n_good_ccd = {}\n'.format(n_good_ccd)
    if apply_subtraction:
        text += 'run {}\n'.format(run)
        text += 'skyscale = {:.1f}\n'.format(skyscale)
    plt.text(1.08, 0.710, text, fontsize=16)

    plt.axis([1.1, -1.1, -1.05, 1.05])
    plt.axis('off')
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    # plt.colorbar(fraction=0.04, pad=0.04)
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()



tt = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11dev/survey-ccds-decam-dr11-trim.fits', columns=['expnum', 'ccd_cuts']))
print('CCD', len(tt))
mask = tt['ccd_cuts']==0
tt = tt[mask]
print('CCD', len(tt))

ss = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/data/dr11dev/sky_scales/skyscales.fits'))
print('skyscales', len(ss))
mask = np.in1d(ss['expnum'], tt['expnum'])
ss = ss[mask]
print('skyscales', len(ss))

# np.random.seed(777)
# idx = np.random.choice(len(ss), size=256)
# expnum_list = np.array(ss['expnum'][idx])
# print(len(expnum_list), 'images to plot')

np.random.seed(666)
runs = np.unique(ss['run'])
ss['random'] = np.random.rand(len(ss))
ss.sort('random')
_, idx_keep = np.unique(ss['run'], return_index=True)
ss = ss[idx_keep]
print('skyscales', len(ss))
expnum_list = np.array(ss['expnum'])
print(len(expnum_list), 'images to plot')

zipped_arg_list = list(zip(expnum_list, repeat(True)))
zipped_arg_list += list(zip(expnum_list, repeat(False)))
with Pool(processes=n_processes) as pool:
    res = pool.starmap(make_plots, zipped_arg_list)

print('Done!!!!!!!!!!!!!!!!!!!!!')
