from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter

from multiprocessing import Pool
import argparse
from pathlib import Path

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'} 
plt.rcParams.update(params)

n_processes = 6
max_exposure = 50
smoothing_scale = 40

parser = argparse.ArgumentParser()
parser.add_argument('n_task')
parser.add_argument('task_id')
args = parser.parse_args()
n_task = int(args.n_task)
task_id = int(args.task_id)

#############################################################################

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

ccdnamenumdict = {'CCD1': 1, 'CCD2': 2, 'CCD3': 3, 'CCD4': 4}
ccdnamenumdict_inv = {aa: bb for bb, aa in ccdnamenumdict.items()}

ccdnum_list = [1, 2, 3, 4]

expnum_blacklist = []

#######################################################################################################################

image_dir = '/global/project/projectdirs/cosmo/staging'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/90prime_ccd_blob_mask'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-90prime-dr9.fits.gz'
output_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates_90prime'

# ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts', 'ccdskycounts', 'plver']
ccd_columns = ['image_hdu', 'expnum', 'ccdname', 'ccdskycounts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
print(len(ccd))

# skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8-90prime.fits')
skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8-90prime-subset.fits')
print(len(skyrun))

mask = skyrun['ok']==True
skyrun = skyrun[mask]
print(len(skyrun))

# # Exclude templates already created -- caution!
# fn_list = glob.glob(os.path.join(output_dir, '*.fits.fz'))
# run_list_done = [int(fn[len(os.path.join(output_dir, 'sky_templates_'))+1:-8]) for fn in fn_list]
# mask = ~np.in1d(skyrun['run'], run_list_done)
# skyrun = skyrun[mask]
# print(len(skyrun), len(run_list_done))

# # Wait to avoid race condition from writing files and checking file status
# time.sleep(60)

# ########################## Exclude z band ##########################
# band = 'z'
# mask = skyrun['filter']!=band
# skyrun = skyrun[mask]
# print(len(skyrun))
# ####################################################################

run_list = np.unique(skyrun['run'])
print(len(run_list))

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
run_list = np.random.choice(run_list, size=len(run_list), replace=False)

# split among the Cori nodes
run_list_split = np.array_split(run_list, n_task)
run_list = run_list_split[task_id]
print('Number of runs in this node:', len(run_list))

#######################################################################################################################

def compute_smooth_sky(run, diagnostic_touch=True):

    skyrun_idx = np.where(skyrun['run']==run)[0]
    band = skyrun['filter'][skyrun_idx[0]]

    print('band: {}, run: {}'.format(band, run))

    #############################################
    # Maybe there's a better way to downselect the exposures?
    if len(skyrun_idx>max_exposure):
        skyrun_idx = skyrun_idx[:max_exposure]
    #############################################

    output_path = os.path.join(output_dir, 'sky_template_{}_{}.fits.fz'.format(band, run))
    residual_path = os.path.join(output_dir, 'sky_residual_{}_{}.fits.fz'.format(band, run))
    if os.path.isfile(output_path):
        print(output_path+' already exists!')
        return None
        
    if diagnostic_touch:
        Path('/global/u2/r/rongpu/temp/sky_template_status/'+os.path.basename(output_path)).touch()
        Path('/global/u2/r/rongpu/temp/sky_template_being_written/'+os.path.basename(output_path)).touch()

    hdul_template = fitsio.FITS(output_path, mode='rw', clobber=True)
    hdul_template.write(data=None) # first HDU is empty

    hdul_residual = fitsio.FITS(residual_path, mode='rw', clobber=True)
    hdul_residual.write(data=None) # first HDU is empty

    for ccdnum in ccdnum_list:

        # print(ccdnum)

        # ####################
        # start = time.time()
        # ####################

        img_list = []
        ccdname = ccdnamenumdict_inv[ccdnum]
        
        for index, skyrun_index in enumerate(skyrun_idx):

            expnum = skyrun['expnum'][skyrun_index]
            
            if expnum in expnum_blacklist:
                continue
            
            # print(ccdnum, ccdname, index, '/', len(skyrun_idx))

            # Load CCD image
            img_fn = os.path.join(image_dir, skyrun['image_filename'][skyrun_index]).strip()
            
            try:
                img = fitsio.read(img_fn, ext=ccdname)
            except:
                print(ccdname+' '+img_fn+' does not exist!')
                continue

            # Get HDU index
            with fitsio.FITS(img_fn) as f:
                hdu_index = f.movnam_ext(ccdname)

            # Load blob mask
            str_loc = str.find(skyrun['image_filename'][skyrun_index].strip(), '.fits')
            img_filename_base = skyrun['image_filename'][skyrun_index].strip()[:str_loc]
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
            
            # Remove median sky
            sky = np.median(img[blob].flatten())
            img = img - sky

            # # Find the entry in survey-ccd
            # ccd_index = np.where((ccd['expnum']==skyrun['expnum'][skyrun_index]) & (ccd['image_hdu']==hdu_index))[0][0]

            # Get the median ccdskycounts of the exposure
            mask = ccd['expnum']==skyrun['expnum'][skyrun_index]
            ccdskycounts_median = np.median(ccd['ccdskycounts'][mask])
            # print('ccdskycounts_median = {:.4f}'.format(ccdskycounts_median))
            
            # Normalize by ccdskycounts
            img = img/ccdskycounts_median
            
            # Apply blob mask
            img[~blob] = np.nan

            img_list.append(img)

            gc.collect()

        if len(img_list)==0:
            print('There is no available {} CCD'.format(ccdname))
            continue

        img_median = np.nanmedian(img_list, axis=0)

        # Fill in NAN values
        mask = ~np.isfinite(img_median)
        # print('number of NAN pixels:', np.sum(mask))
        img_median[mask] = 0

        img_median1 = img_median.copy()

        # 3-sigma clipping
        sky_nmad = nmad(img_median[np.isfinite(img_median)]) # sky level
        # mask = (img_median<-3*sky_nmad) | (img_median>3*sky_nmad)
        # img_median1[mask] = 0
        mask = (img_median<-3*sky_nmad)
        img_median1[mask] = -3*sky_nmad
        mask = (img_median>3*sky_nmad)
        img_median1[mask] = 3*sky_nmad

        # trim edges
        trim_size = 10            
        img_median1 = img_median1[trim_size:(img_median1.shape[0]-trim_size), trim_size:(img_median1.shape[1]-trim_size)]

        # downsize the image to speed up gaussian filter
        binsize = 1 # no downsizing
        # img_median1 = np.nanmean(np.nanmean(img_median1.reshape((img_median1.shape[0]//binsize, binsize, img_median1.shape[1]//binsize,-1)), axis=3), axis=1)
        x_small_grid = trim_size + binsize/2 + binsize*np.arange(img_median1.shape[1])
        y_small_grid = trim_size + binsize/2 + binsize*np.arange(img_median1.shape[0])

        # Gaussian filtering
        img_median1_smooth = gaussian_filter(img_median1, smoothing_scale, mode='reflect')

        # Convert the downsized smooth image to full size
        interp_func = interp2d(x_small_grid, y_small_grid, img_median1_smooth, kind='linear')
        x_grid, y_grid = np.arange(img_median.shape[1]), np.arange(img_median.shape[0])
        img_median_smooth = interp_func(x_grid, y_grid).reshape(img_median.shape)

        ################################ Save sky template ###############################

        hdul_template.write(data=img_median_smooth, extname=ccdname, compress='rice')
        img_median -= img_median_smooth
        hdul_residual.write(data=img_median, extname=ccdname, compress='rice')

        # ##################
        # end = time.time()
        # print('Took {:.1f} seconds'.format(end - start))
        # ##################

    hdul_template.close()
    hdul_residual.close()
    
    if diagnostic_touch:
        os.remove('/global/u2/r/rongpu/temp/sky_template_being_written/'+os.path.basename(output_path))

def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(compute_smooth_sky, run_list)

    print('All done!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

