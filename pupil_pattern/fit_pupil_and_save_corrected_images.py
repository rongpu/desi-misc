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

n_processes = 32

cp_image_dir = '/global/project/projectdirs/cosmo/staging'
lg_image_dir = '/global/cfs/cdirs/cosmo/staging/decam/CP-LG9'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'

template_dir = '/global/cscratch1/sd/rongpu/dr9dev/pupil_pattern/pupil_templates/'
output_dir = '/global/cscratch1/sd/rongpu/dr9dev/pupil_pattern/pupil_corrected_images/'

pix_size = 0.262/3600

# image_path_list = glob.glob(os.path.join(cp_image_dir, '*ooi*.fits.fz'))
ccd = Table.read('/global/cfs/cdirs/desi/users/rongpu/dr9/pupil_pattern/survey-ccds-decam-dr9-lg-pupilrun.fits')
# mask = ccd['blacklist']==False
# ccd = ccd[mask]

image_fn_lg_list = np.unique(ccd['image_filename_lg'])

def fit_pupil(image_fn_lg):
    
    mask = ccd['image_filename_lg']==image_fn_lg
    run = ccd['run'][mask][0]
    band = ccd['filter'][mask][0]
    cp_image_filename = ccd['image_filename'][mask][0]

    print('image: {}, band: {}, run: {}'.format(image_fn_lg, band, run))

    pupil_path = os.path.join(template_dir, 'pupil_smooth_{}_{}.fits.fz'.format(band, run))
    lg_image_path = os.path.join(lg_image_dir, image_fn_lg)
    cp_img_fn = os.path.join(cp_image_dir, cp_image_filename).strip()
    ood_fn = cp_img_fn.replace('_ooi_', '_ood_')
    output_path = os.path.join(output_dir, image_fn_lg)

    if os.path.isfile(output_path):
        print(output_path+' already exists!')
        return None

    ######################################## Get the scaling factor ########################################

    median_image_list = []
    median_pupil_list = []

    for ii, ccdnum in enumerate(ccdnum_list):

        ccdname = ccdnamenumdict_inv[ccdnum]

        if ccdname=='S7':
            continue

        ccd_index = np.where((ccd['image_filename_lg']==image_fn_lg) & (ccd['ccdname']==ccdname))[0]
        if len(ccd_index)==1:
            ccd_index = ccd_index[0]
        elif len(ccd_index)>1:
            raise ValueError
        elif len(ccd_index)==0:
            continue
        
        expnum = ccd['expnum'][ccd_index]

        bestfit_angle = ccd['gradient_angle'][ccd_index]
        bestfit_slope = ccd['gradient_slope'][ccd_index]
        bestfit_intercept = ccd['gradient_intercept'][ccd_index]
        bestfit_unit_x, bestfit_unit_y = np.cos(bestfit_angle/180.*np.pi), np.sin(bestfit_angle/180.*np.pi)
                    
        # print(ccdnum, ccdname, index, '/', len(ccd_idx))

        # Load CCD image
        lg_image_path = os.path.join(lg_image_dir, ccd['image_filename_lg'][ccd_index]).strip()
        
        try:
            pupil = fitsio.read(pupil_path, ext=ccdname)
        except:
            print(ccdname+' '+pupil_path+' does not exist!')
            continue

        try:
            img = fitsio.read(lg_image_path, ext=ccdname)
        except:
            print(ccdname+' '+lg_image_path+' does not exist!')
            continue

        # Apply gradient correction
        ra, dec = ccd_ra[ii], ccd_dec[ii]
        pix_y_grid, pix_x_grid = np.meshgrid(np.arange(img_shape[1]), np.arange(img_shape[0]))
        ra_pix = ra + (pix_x_grid-img_shape[0]/2) * pix_size
        dec_pix = dec - (pix_y_grid-img_shape[1]/2) * pix_size
        x, y = ra_pix, dec_pix
        dotprod_list = x * bestfit_unit_x + y * bestfit_unit_y
        img -= bestfit_slope * dotprod_list + bestfit_intercept

        try:
            ood = fitsio.read(ood_fn, ext=ccdname)
        except:
            print(ccdname+' '+ood_fn+' does not exist!')
            continue

        # Get HDU index
        with fitsio.FITS(cp_img_fn) as f:
            hdu_index = f.movnam_ext(ccdname)

        # Load blob mask
        str_loc = str.find(ccd['image_filename'][ccd_index].strip(), '.fits')
        img_filename_base = ccd['image_filename'][ccd_index].strip()[:str_loc]
        blob_path = os.path.join(blob_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
        try:
            blob_data = np.load(blob_path)
        except:
            #################
            # DO SOMETHING HERE?
            #################
            print(blob_path+' does not exist!')
            continue
        try:
            blob = blob_data['hdu'+str(hdu_index).zfill(2)]
        except:
            print(blob_path+' hdu'+str(hdu_index)+' does not exist!')
            continue

        # # Subtract sky level
        # img -= bestfit_intercept

        # Apply blob and ood mask
        img_mask = (blob==True) & (ood==0)
        img[~img_mask] = np.nan

        if np.sum(np.isfinite(img)) < 0.2 * np.prod(img_shape):
            print('not enough pixels for sky estimate: only {:.1f}% avaialbe'.format(np.sum(np.isfinite(img))/np.prod(img_shape)*100))
            continue

        median_image_list.append(np.median(img[np.isfinite(img)]))
        median_pupil_list.append(np.median(pupil))

    pupil_slope, pupil_intercept, r_value, p_value, std_err = stats.linregress(median_pupil_list, median_image_list)
    print('slope: {}, intercept: {}'.format(pupil_slope, pupil_intercept))

    ######################################## Apply pupil correction ########################################

    hdul = fitsio.FITS(output_path, mode='rw', clobber=True)
    hdul.write(data=None) # first HDU is empty

    for ii, ccdnum in enumerate(ccdnum_list):

        ccdname = ccdnamenumdict_inv[ccdnum]

        ccd_index = np.where((ccd['image_filename_lg']==image_fn_lg) & (ccd['ccdname']==ccdname))[0]
        if len(ccd_index)==1:
            ccd_index = ccd_index[0]
        elif len(ccd_index)>1:
            raise ValueError
        elif len(ccd_index)==0:
            continue
        
        expnum = ccd['expnum'][ccd_index]

        bestfit_angle = ccd['gradient_angle'][ccd_index]
        bestfit_slope = ccd['gradient_slope'][ccd_index]
        bestfit_intercept = ccd['gradient_intercept'][ccd_index]
        bestfit_unit_x, bestfit_unit_y = np.cos(bestfit_angle/180.*np.pi), np.sin(bestfit_angle/180.*np.pi)
                    
        # print(ccdnum, ccdname, index, '/', len(ccd_idx))

        # Load CCD image
        lg_image_path = os.path.join(lg_image_dir, ccd['image_filename_lg'][ccd_index]).strip()
        
        try:
            pupil = fitsio.read(pupil_path, ext=ccdname)
        except:
            print(ccdname+' '+pupil_path+' does not exist!')
            continue

        try:
            img = fitsio.read(lg_image_path, ext=ccdname)
        except:
            print(ccdname+' '+lg_image_path+' does not exist!')
            continue

        # Apply gradient correction
        ra, dec = ccd_ra[ii], ccd_dec[ii]
        pix_y_grid, pix_x_grid = np.meshgrid(np.arange(img_shape[1]), np.arange(img_shape[0]))
        ra_pix = ra + (pix_x_grid-img_shape[0]/2) * pix_size
        dec_pix = dec - (pix_y_grid-img_shape[1]/2) * pix_size
        x, y = ra_pix, dec_pix
        dotprod_list = x * bestfit_unit_x + y * bestfit_unit_y
        img -= bestfit_slope * dotprod_list

        img -= pupil_slope * pupil

        hdul.write(data=img, extname=ccdname, compress='rice')

    hdul.close()
    
def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(fit_pupil, image_fn_lg_list)
    
    # make_plots('c4d_181102_051559_ooi_z_lg9.fits.fz')

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

