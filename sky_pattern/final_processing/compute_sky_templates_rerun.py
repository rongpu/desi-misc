# Rerun with bad exposures excluded

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

n_processes = 12

parser = argparse.ArgumentParser()
parser.add_argument('n_task')
parser.add_argument('task_id')
args = parser.parse_args()
n_task = int(args.n_task)
task_id = int(args.task_id)

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
                  'N29': 60, 'N30': 61, 'N31': 62,
                  }
ccdnamenumdict_inv = {aa: bb for bb, aa in ccdnamenumdict.items()}

ccdnum_list = [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
       18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
       35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
       52, 53, 54, 55, 56, 57, 58, 59, 60, 62]

# Shape of the DECam CP image
img_shape = (4094, 2046)

max_exposure = 50

expnum_blacklist = [243224, 243233, 243250, 243261, 247512, 247519, 247524, 247535,
       247536, 247543, 247544, 247546, 251698, 261222, 261236, 261245,
       261264, 261289, 263153, 263682, 269550, 269553, 269556, 269573,
       269575, 269583, 269584, 269662, 269669, 269719, 270231, 270260,
       276448, 276450, 276454, 449966, 463823, 569657, 600963, 600966,
       611450, 690302, 718586, 718600, 718608, 718626, 718636, 720061,
       720100, 754068, 754083, 768677, 803349, 803356, 807350, 808254,
       808338, 808339, 808547, 808560, 808570, 808650, 863383]

#######################################################################################################################

# blob_dir = '/global/cscratch1/sd/rongpu/fringe/decam_ccd_blob_mask'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'

image_dir = '/global/project/projectdirs/cosmo/staging'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
output_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates_rerun'

# ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts', 'ccdskycounts', 'plver']
ccd_columns = ['image_hdu', 'expnum', 'ccdname', 'ccdskycounts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
print(len(ccd))

skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8.fits')
# skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8_less.fits')
print(len(skyrun))

mask = skyrun['ok']==True
skyrun = skyrun[mask]
print(len(skyrun))

# # Exclude templates already created
# fn_list = glob.glob(os.path.join(output_dir, '*.fits.fz'))
# run_list_done = [int(fn[len(os.path.join(output_dir, 'sky_templates_'))+1:-8]) for fn in fn_list]
# mask = ~np.in1d(skyrun['run'], run_list_done)
# skyrun = skyrun[mask]
# print(len(skyrun), len(run_list_done))

# ########################## Exclude z band ##########################
# band = 'z'
# mask = skyrun['filter']!=band
# skyrun = skyrun[mask]
# print(len(skyrun))
# ####################################################################

# run_list = np.unique(skyrun['run'])
# print(len(run_list))

run_list = np.array([16, 20, 31, 34, 126, 249, 342, 345, 369, 457, 459, 791, 839, 855, 879])

# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
run_list = np.random.choice(run_list, size=len(run_list), replace=False)

# split among the Cori nodes
run_list_split = np.array_split(run_list, n_task)
run_list = run_list_split[task_id]
print('Number of runs in this node:', len(run_list))

# Wait to avoid race condition from writing files and checking file status
time.sleep(60)

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
    if os.path.isfile(output_path):
        print(output_path+' already exists!')
        return None
        
    if diagnostic_touch:
        Path('/global/u2/r/rongpu/temp/sky_template_status/'+os.path.basename(output_path)).touch()
        Path('/global/u2/r/rongpu/temp/sky_template_being_written/'+os.path.basename(output_path)).touch()

    hdul_w = fitsio.FITS(output_path, mode='rw', clobber=True)
    hdul_w.write(data=None) # first HDU is empty

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
            
            if ccdname=='S7':
                # Only keep the good half of the S7
                half = img_shape[1] // 2
                img = img[:, :half]
                blob = blob[:, :half]

            # Remove median sky
            sky = np.median(img[blob].flatten())
            img = img - sky

            # # Find the entry in survey-ccd
            # ccd_index = np.where((ccd['expnum']==skyrun['expnum'][skyrun_index]) & (ccd['image_hdu']==hdu_index))[0][0]

            # Get the median ccdskycounts of the exposure
            # mask = (ccd['expnum']==skyrun['expnum'][skyrun_index]) & (ccd['ccdname']!='S7') & (ccd['ccdname']!='S7 ') # too slow!
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

        if ccdname=='S7':

            # trim three of edges
            trim_size = 20
            trim_size_top = 1
            img_median1 = img_median1[trim_size:(img_median1.shape[0]-trim_size), trim_size_top:(img_median1.shape[1]-trim_size)]

            # downsize the image to speed up gaussian filter
            binsize = 2
            img_median1 = np.nanmean(np.nanmean(img_median1.reshape((img_median1.shape[0]//binsize, binsize, img_median1.shape[1]//binsize,-1)), axis=3), axis=1)
            x_small_grid = trim_size_top + binsize/2 + binsize*np.arange(img_median1.shape[1])
            y_small_grid = trim_size + binsize/2 + binsize*np.arange(img_median1.shape[0])

        elif (band=='r' and skyrun['expnum'][skyrun_index]<298251):

            ################### r band edge glow exposures ###################

            # r-band edge glow
            # trim edges
            trim_size = 10
            extra_trim_size = 50
            img_median1 = img_median1[trim_size:(img_median1.shape[0]-trim_size), trim_size:(img_median1.shape[1]-trim_size)]

            # The edge glow only appears on one edge and are at opposite edges for top CCDs and bottom CCDs
            if ccdnum<=31:
                img_median1 = img_median1[:-extra_trim_size]
            else:
                img_median1 = img_median1[extra_trim_size:]

            # downsize the image to speed up gaussian filter
            binsize = 2
            img_median1 = np.nanmean(np.nanmean(img_median1.reshape((img_median1.shape[0]//binsize, binsize, img_median1.shape[1]//binsize,-1)), axis=3), axis=1)

            x_small_grid = trim_size + binsize/2 + binsize*np.arange(img_median1.shape[1])
            if ccdnum<=31:
                y_small_grid = trim_size + binsize/2 + binsize*np.arange(img_median1.shape[0])
            else:
                y_small_grid = trim_size + extra_trim_size + binsize/2 + binsize*np.arange(img_median1.shape[0])

        else:

            ################### Normal exposures ####################

            # trim edges
            trim_size = 10            
            img_median1 = img_median1[trim_size:(img_median1.shape[0]-trim_size), trim_size:(img_median1.shape[1]-trim_size)]

            # downsize the image to speed up gaussian filter
            binsize = 2
            img_median1 = np.nanmean(np.nanmean(img_median1.reshape((img_median1.shape[0]//binsize, binsize, img_median1.shape[1]//binsize,-1)), axis=3), axis=1)
            x_small_grid = trim_size + binsize/2 + binsize*np.arange(img_median1.shape[1])
            y_small_grid = trim_size + binsize/2 + binsize*np.arange(img_median1.shape[0])

        # Gaussian filtering
        img_median1_smooth = gaussian_filter(img_median1, 60, mode='reflect')

        # Convert the downsized smooth image to full size
        interp_func = interp2d(x_small_grid, y_small_grid, img_median1_smooth, kind='linear')
        x_grid, y_grid = np.arange(img_median.shape[1]), np.arange(img_median.shape[0])
        img_median_smooth = interp_func(x_grid, y_grid).reshape(img_median.shape)

        if ccdname=='S7':
            # Add back the other half
            tmp = np.zeros(img_shape)
            half = img_shape[1] // 2
            tmp[:, :half] = img_median_smooth
            img_median_smooth = tmp


        ################################ Save sky template ###############################

        hdul_w.write(data=img_median_smooth, extname=ccdname, compress='rice')

        # ##################
        # end = time.time()
        # print('Took {:.1f} seconds'.format(end - start))
        # ##################

    hdul_w.close()
    
    if diagnostic_touch:
        os.remove('/global/u2/r/rongpu/temp/sky_template_being_written/'+os.path.basename(output_path))

def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(compute_smooth_sky, run_list)

    print('All done!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

