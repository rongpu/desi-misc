# Rerun template fitting for a subset of the runs

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
# matplotlib.use('Agg')
# Use the same backend as Jupyter notebook
matplotlib.use('module://ipykernel.pylab.backend_inline')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
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

ccd_dec = [0.90299039, 0.90274404, 0.90285652, 0.73894001, 0.73933177, 0.73919444,
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

run_list = [436, 437, 438, 439, 445, 446, 447, 455, 456, 469, 471]

overwrite = True

################################################################################

n_processes = 32

parser = argparse.ArgumentParser()
parser.add_argument('n_task')
parser.add_argument('task_id')
args = parser.parse_args()
n_task = int(args.n_task)
task_id = int(args.task_id)

image_dir = '/global/project/projectdirs/cosmo/staging'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# template_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates/'
template_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates_final_edge_glow'
########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
skyscale_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales/'

skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8.fits')
print('skyrun', len(skyrun))

skyrun_all = skyrun.copy()

# sky_path_list = glob.glob(os.path.join(template_dir, '*.fits.fz'))

# ####################################################################################
# # The file should be at least 5 hours old to ensure it's not being written
# for sky_path in sky_path_list:
#     time_modified = os.path.getmtime(sky_path)
#     if (time.time() - time_modified)/3600 < 5:
#         sky_path_list.remove(sky_path)
# print('sky_path_list', len(sky_path_list))
# ####################################################################################

# #################################### Exclude z band ####################################
# mask = np.array(['_z_' in sky_path for sky_path in sky_path_list])
# sky_path_list = np.array(sky_path_list)[~mask]
# ########################################################################################

# run_list = np.array([int(fn[len(os.path.join(template_dir, 'sky_templates_'))+1:-8]) for fn in sky_path_list])
# print('run_list', len(run_list))

mask = np.in1d(skyrun['run'], run_list)
skyrun = skyrun[mask]
print('skyrun', len(skyrun))

# # Remove completed runs from list
# expnum_status = Table.read('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/fitting_status.fits')
# expnum_status = expnum_status[expnum_status['done']==False]
# mask = np.in1d(skyrun['expnum'], expnum_status['expnum'])
# skyrun = skyrun[mask]
# print('skyrun', len(skyrun))

expnum_list = skyrun['expnum'].copy()

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)

# split among the Cori nodes
expnum_list_split = np.array_split(expnum_list, n_task)
expnum_list = expnum_list_split[task_id]
print('Number of exposures in this node:', len(expnum_list))

ccd_columns = ['image_hdu', 'expnum', 'ccdname', 'ccdskycounts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))

plot_dir = '/global/cfs/cdirs/desi/www/users/rongpu/plots/dr9dev/sky_pattern/sky_templates_v2/fit_scale'

binsize = 2
pix_size = 0.262/3600*binsize

plot_q = True
plots_per_run = 3

image_vrange = {'g':5, 'r':6, 'z':30}

def template_fitting(expnum, diagnostic_touch=True):
    
    # # The file should be at least 5 hours old to ensure it's not being written
    # time_modified = os.path.getmtime(sky_path)
    # if (time.time() - time_modified)/3600 < 5:
    #     # continue
    #     return None

    # run = int(sky_path[len(os.path.join(template_dir, 'sky_templates_'))+1:-8])

    # Get run info
    skyrun_index = np.where(skyrun['expnum']==expnum)[0][0]
    band = skyrun['filter'][skyrun_index]
    run = skyrun['run'][skyrun_index]
    
    sky_path = os.path.join(template_dir, 'sky_template_{}_{}.fits.fz'.format(band, run))

    vrange = image_vrange[band]

    mask = skyrun['run']==run
    skyrun_idx = np.where(mask)[0]
    # print('\nrun {}, {} exposures'.format(run, len(skyrun_idx)))

    mask = skyrun_all['run']==run
    expnum_list_plot = skyrun_all['expnum']
    
    np.random.seed(123+run)
    expnum_list_plot = np.random.choice(expnum_list_plot, size=plots_per_run, replace=False)

    ####################
    start = time.time()
    ####################

    image_filename = skyrun['image_filename'][skyrun_index].strip()
    image_path = os.path.join(image_dir, image_filename)
    ood_path = image_path.replace('_ooi_', '_ood_')

    blob_path = os.path.join(blob_dir, 'blob_mask', image_filename.replace('.fits.fz', '-blobmask.npz'))
    skyscale_path = os.path.join(skyscale_dir, image_filename.replace('.fits.fz', '-skyscale.txt'))

    if not os.path.exists(os.path.dirname(skyscale_path)):
        try: # in case another process is also creating the directory
            os.makedirs(os.path.dirname(skyscale_path))
        except:
            pass

    if (overwrite==False) and os.path.isfile(skyscale_path):
        # print(skyscale_path+' already exists!!')
        return None

    result = Table(names=('image_hdu', 'ccdname', 'ccdskyscale', 'medianskyscale'), dtype=('i4', 'S3', 'f4', 'f4'))
    result['ccdskyscale'].format = '.5f'
    result['medianskyscale'].format = '.5f'

    # try:
    #     blob_data = np.load(blob_path)
    # except:
    #     print(blob_path+' does not exist!')
    #     result.write(skyscale_path, format='ascii.commented_header', overwrite=True)
    #     continue

    if os.stat(blob_path).st_size == 0:
        print(blob_path+' is empty!')
        result.write(skyscale_path, format='ascii.commented_header', overwrite=True)
        return None
    else:
        blob_data = np.load(blob_path)

    Path(skyscale_path).touch()
    if diagnostic_touch:
        Path('/global/u2/r/rongpu/temp/sky_scale_being_written/expnum_{}'.format(expnum)).touch()

    if plot_q and (expnum in expnum_list_plot):

        plot_path = os.path.join(plot_dir, band, '{}_{}_image_{}_fitscale.png'.format(band, run, expnum))

        if (overwrite==False) and os.path.isfile(plot_path):
            print(plot_path, 'already exists!!! overwrite')
            # continue

        if not os.path.exists(os.path.dirname(plot_path)):
            try: # in case another process is also creating the directory
                os.makedirs(os.path.dirname(plot_path))
            except:
                pass

        scale_min, scale_max = np.inf, -np.inf
        scale_list = []

        plt.figure(figsize=(13.7*2, 13.075*2))

    for ii, ccdnum in enumerate(ccdnum_list):

        ccdname = ccdnamenumdict_inv[ccdnum]
        # print(ii, ccdname)

        if ccdname=='S7':
            continue

        try:
            img = fits.getdata(image_path, extname=ccdname)
            ood = fits.getdata(ood_path, extname=ccdname)
        except:
            print(ccdname+' does not exist in image!')
            continue

        try:
            sky = fits.getdata(sky_path, extname=ccdname)
        except:
            print(ccdname+' does not exist in template!')
            continue

        if np.all(sky==0):
            print(ccdname+' template is all zeros!')
            continue

        # Find the entry in survey-ccd
        if len(ccdname)==3:
            ccdname_space_filled = ccdname
        else:
            ccdname_space_filled = ccdname+' '
        ccd_index = np.where((ccd['expnum']==expnum) & (ccd['ccdname']==ccdname_space_filled))[0][0]

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
            print(blob_path+' hdu'+str(hdu_index)+' does not exist!')
            continue

        if ccdname=='S7':
            # Only keep the good half of the S7
            half = img_shape[1] // 2
            img = img[:, :half]
            ood = ood[:, :half]
            sky = sky[:, :half]
            blob = blob[:, :half]

        # # Remove median sky
        # sky = np.median(img[blob].flatten())
        # img = img - sky

        # naive sky estimation
        mask = (img<np.percentile(img.flatten(), 95))
        median_sky = np.median(img[mask].flatten())
        img = img - median_sky

        # Apply blob mask
        img1 = img.copy()
        img1_mask = (blob==True) & (ood==0)
        img1[~img1_mask] = np.nan
        img1_mask = np.isfinite(img1)

        if np.sum(~img1_mask)/np.prod(img1_mask.shape)>0.8:
            print('{} {:.3f}% pixels are masked; skip'.format(ccdname, np.sum(~img1_mask)/np.prod(img1_mask.shape)*100))
            continue
        # elif np.sum(~img1_mask)/np.prod(img1_mask.shape)>0.3:
        #     print('{} {:.3f}% pixels are masked'.format(ccdname, np.sum(~img1_mask)/np.prod(img1_mask.shape)*100))

        # 3-sigma clipping
        img1_nmad = nmad(img1[img1_mask]) # sky level
        # mask = (img1<-3*img1_nmad) | (img1>3*img1_nmad)
        # img1[mask] = 0
        ##############################################################
        mask = (img1<-3*img1_nmad) | (img1>3*img1_nmad)
        if np.sum(mask)/np.sum(img1_mask)>0.03:
            print('{} {:.3f}% pixels are clipped'.format(ccdname, np.sum(mask)/np.sum(img1_mask)*100))
        ##############################################################
        mask = (img1<-3*img1_nmad)
        img1[mask] = -3*img1_nmad
        mask = (img1>3*img1_nmad)
        img1[mask] = 3*img1_nmad

        img1_flat = img1[img1_mask].flatten()
        sky_flat = sky[img1_mask].flatten()
        slope, intercept, r_value, p_value, std_err = stats.linregress(sky_flat, img1_flat)
        result.add_row((hdu_index, ccdname, slope, 0))

        if plot_q and (expnum in expnum_list_plot):

            if (slope < scale_min) and (ccdname!='S7'):
                scale_min = slope
                scale_min_ccdname = ccdname
            if (slope > scale_max) and (ccdname!='S7'):
                scale_max = slope
                scale_max_ccdname = ccdname

            # if slope>2 or slope<0:
            #     print('{} slope, intercept = {:.4f}, {:.4f}'.format(ccdname, slope, intercept))

            # Apply sky pattern correction
            img = img - sky * slope

            if ccdname=='S7':
                # Add back the other half
                tmp = np.zeros(img_shape)
                half = img_shape[1] // 2
                tmp[:, :half] = img
                img = tmp

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

            fig = plt.imshow(img.T, cmap='seismic', vmin=-vrange, vmax=vrange, 
                       extent=(ra-ysize*pix_size/2, ra+ysize*pix_size/2, dec-xsize*pix_size/2, dec+xsize*pix_size/2))

    if len(result)==0:
        n_ccd = 0
    else:
        mask = result['ccdname']!='S7'
        n_ccd = np.sum(mask)
    if n_ccd>0:
        medianskyscale = np.median(result['ccdskyscale'][mask])
        result['medianskyscale'] = medianskyscale
    else:
        print('No CCD available!')
    result.write(skyscale_path, format='ascii.commented_header', overwrite=True)

    if diagnostic_touch:
        os.remove('/global/u2/r/rongpu/temp/sky_scale_being_written/expnum_{}'.format(expnum))

    if plot_q and (expnum in expnum_list_plot) and (len(result)>0):

        print('making plots')

        text = 'run {}, {} band\n'.format(run, band)
        text += 'expnum = {}\n'.format(expnum)
        text += 'scale = per-CCD fit\n'
        text += 'median scale = {:.1f}\n'.format(result['medianskyscale'][0])
        text += 'min scale = {:.1f} ({})\n'.format(scale_min, scale_min_ccdname)
        text += 'max scale = {:.1f} ({})\n'.format(scale_max, scale_max_ccdname)
        plt.text(1.08, 0.710, text, fontsize=16)

        plt.axis([1.1, -1.1, -1.05, 1.05])
        plt.axis('off')
        fig.axes.get_xaxis().set_visible(False)
        fig.axes.get_yaxis().set_visible(False)
        # plt.colorbar(fraction=0.04, pad=0.04)
        plt.tight_layout()
        plt.savefig(plot_path)
        plt.close()

    ##################
    end = time.time()
    print('Took {:.1f} seconds'.format(end - start))
    ##################

def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(template_fitting, expnum_list)

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()
