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
from scipy import stats
from multiprocessing import Pool
import argparse

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'} 
plt.rcParams.update(params)

n_processes = 32
overwrite = False
plot_q = True
binsize = 2
pix_size = 0.262/3600*binsize
######################
# plots_per_run = 3
plots_per_run = np.inf
######################


################################################################################

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

ccdnamenumdict = {'CCD1': 1, 'CCD2': 2, 'CCD3': 3, 'CCD4': 4}
ccdnamenumdict_inv = {aa: bb for bb, aa in ccdnamenumdict.items()}

ccdnum_list = [1, 2, 3, 4]
ccd_ra = [-0.1554, -0.1554, 0.1554, 0.1554]
ccd_dec = [-0.1554, 0.1554, -0.1554, 0.1554]

################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('n_task')
parser.add_argument('task_id')
args = parser.parse_args()
n_task = int(args.n_task)
task_id = int(args.task_id)

image_dir = '/global/project/projectdirs/cosmo/staging'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-mosaic-dr9.fits.gz'
########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
template_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates_mosaic'
########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!########################!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
skyscale_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales_mosaic/'

# skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8-mosaic.fits')
skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8-mosaic-subset.fits')
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

image_vrange = {'g':0.25, 'r':0.25, 'z':0.4}

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

    if np.isfinite(plots_per_run):
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

        plt.figure(figsize=(8, 8))
        # for 90prime:
        # plt.figure(figsize=(11, 11))

    for ii, ccdnum in enumerate(ccdnum_list):

        ccdname = ccdnamenumdict_inv[ccdnum]
        # print(ii, ccdname)

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

            if (slope < scale_min):
                scale_min = slope
                scale_min_ccdname = ccdname
            if (slope > scale_max):
                scale_max = slope
                scale_max_ccdname = ccdname

            # if slope>2 or slope<0:
            #     print('{} slope, intercept = {:.4f}, {:.4f}'.format(ccdname, slope, intercept))

            # Apply sky pattern correction
            img = img - sky * slope

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
                       extent=(ra-ysize*pix_size/2, ra+ysize*pix_size/2, dec+xsize*pix_size/2, dec-xsize*pix_size/2))

    n_ccd = len(result)
    if n_ccd>0:
        medianskyscale = np.median(result['ccdskyscale'][mask])
        result['medianskyscale'] = medianskyscale
    else:
        print('No CCD available!')
    result.write(skyscale_path, format='ascii.commented_header', overwrite=True)

    if diagnostic_touch:
        os.remove('/global/u2/r/rongpu/temp/sky_scale_being_written/expnum_{}'.format(expnum))

    if plot_q and (expnum in expnum_list_plot) and (len(result)>0):

        # print('making plots')

        text = 'run {}, {} band, expnum {}, scale = per-CCD fit\n'.format(run, band, expnum)
        text += 'median scale = {:.1f},  min scale = {:.1f} ({}),  max scale = {:.1f} ({})'.format(result['medianskyscale'][0], scale_min, scale_min_ccdname, scale_max, scale_max_ccdname)
        plt.title(text)
        plt.axis([0.3, -0.3, -0.3, 0.3])
        # for 90prime:
        # plt.axis([0.55, -0.55, -0.55, 0.55])
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
