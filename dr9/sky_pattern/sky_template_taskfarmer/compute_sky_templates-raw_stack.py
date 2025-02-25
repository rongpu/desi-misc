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

smoothing_scale = 120 # in pixels

expnum_blacklist = [243224, 243233, 243250, 243261, 247512, 247519, 247524, 247535,
       247536, 247543, 247544, 247546, 251698, 261222, 261236, 261245,
       261264, 261289, 263153, 263682, 269550, 269553, 269556, 269573,
       269575, 269583, 269584, 269662, 269669, 269719, 270231, 270260,
       276448, 276450, 276454, 449966, 463823, 569657, 600963, 600966,
       611450, 690302, 718586, 718600, 718608, 718626, 718636, 720061,
       720100, 754068, 754083, 768677, 803349, 803356, 807350, 808254,
       808338, 808339, 808547, 808560, 808570, 808650, 863383]

halfed_n10_run_list = [376, 377, 378, 384, 385, 386, 798, 799, 800, 806, 807, 808, 1197, 1198, 1199, 1200, 1206, 1207]

#######################################################################################################################

# blob_dir = '/global/cscratch1/sd/rongpu/fringe/decam_ccd_blob_mask'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'

image_dir = '/global/project/projectdirs/cosmo/staging'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
output_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates'

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

# Exclude templates already created
fn_list = glob.glob(os.path.join(output_dir, '*.fits.fz'))
run_list_done = [int(fn[len(os.path.join(output_dir, 'sky_templates_'))+1:-8]) for fn in fn_list]
mask = ~np.in1d(skyrun['run'], run_list_done)
skyrun = skyrun[mask]
print(len(skyrun), len(run_list_done))

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

# Wait to avoid race condition from writing files and checking file status
time.sleep(60)

#######################################################################################################################

def compute_raw_sky(run, diagnostic_touch=True):

    skyrun_idx = np.where(skyrun['run']==run)[0]
    band = skyrun['filter'][skyrun_idx[0]]

    print('band: {}, run: {}'.format(band, run))

    #############################################
    # Maybe there's a better way to downselect the exposures?
    if len(skyrun_idx>max_exposure):
        skyrun_idx = skyrun_idx[:max_exposure]
    #############################################

    raw_path = os.path.join(output_dir, 'sky_raw_{}_{}.fits.fz'.format(band, run))

    if os.path.isfile(raw_path):
        print(raw_path+' already exists!')
        return None
        
    if diagnostic_touch:
        Path('/global/u2/r/rongpu/temp/sky_template_status/'+os.path.basename(raw_path)).touch()
        Path('/global/u2/r/rongpu/temp/sky_template_being_written/'+os.path.basename(raw_path)).touch()


    hdul_raw_stacked = fitsio.FITS(raw_path, mode='rw', clobber=True)
    hdul_raw_stacked.write(data=None) # first HDU is empty

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
            ood_fn = img_fn.replace('_ooi_', '_ood_')
            
            try:
                img = fitsio.read(img_fn, ext=ccdname)
            except:
                print(ccdname+' '+img_fn+' does not exist!')
                continue

            try:
                ood = fitsio.read(ood_fn, ext=ccdname)
            except:
                print(ccdname+' '+ood_fn+' does not exist!')
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
            
            if (ccdname=='S7') or ((run in halfed_n10_run_list) and (ccdname=='N10')):
                # Only keep the good half of the S7
                half = img_shape[1] // 2
                img = img[:, :half]
                ood = ood[:, :half]
                blob = blob[:, :half]

            # Apply blob and ood mask
            img_mask = (blob==True) & (ood==0)

            # Remove median sky
            sky = np.median(img[img_mask].flatten())
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
            
            img[~img_mask] = np.nan

            img_list.append(img)

            gc.collect()

        if len(img_list)==0:
            print('There is no available {} CCD'.format(ccdname))
            continue

        if len(img_list)<15:
            print('There are only {} images available for {} CCD; skip'.format(len(img_list), ccdname))
            continue

        img_median = np.nanmedian(img_list, axis=0)

        if (ccdname=='S7') or ((run in halfed_n10_run_list) and (ccdname=='N10')):
            # Add back the other half
            tmp = np.zeros(img_shape)
            half = img_shape[1] // 2
            tmp[:, :half] = img_median
            img_median = tmp

        # Fill in NAN values
        mask = ~np.isfinite(img_median)
        # print('number of NAN pixels:', np.sum(mask))
        img_median[mask] = 0

        ################################ Save sky template ###############################

        hdul_raw_stacked.write(data=img_median, extname=ccdname, compress='rice')

        # ##################
        # end = time.time()
        # print('Took {:.1f} seconds'.format(end - start))
        # ##################

    hdul_raw_stacked.close()
    
    if diagnostic_touch:
        os.remove('/global/u2/r/rongpu/temp/sky_template_being_written/'+os.path.basename(raw_path))

def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(compute_raw_sky, run_list)

    print('All done!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

