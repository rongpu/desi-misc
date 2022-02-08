## EXPERIMENTAL! NOT USED!!! ##

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
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

n_processes = 32
diagnostic_touch = True


# run = 0

parser = argparse.ArgumentParser()
parser.add_argument('run')
args = parser.parse_args()
run = int(args.run)

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

ccdnamenumdict = {'S1': 25, 'S2': 26, 'S3': 27, 'S4': 28,
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

ccdnum_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,
               18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
               35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
               52, 53, 54, 55, 56, 57, 58, 59, 60, 62]

# Shape of the DECam CP image
img_shape = (4094, 2046)
binsize = 4  # 4x4 downsizing
# img_shape_downsized = (1024, 512)

max_n_exposure = 250
min_n_exposure = 230
# max_n_exposure = 30
# min_n_exposure = 20

expnum_blacklist = [243224, 243233, 243250, 243261, 247512, 247519, 247524, 247535,
                    247536, 247543, 247544, 247546, 251698, 261222, 261236, 261245,
                    261264, 261289, 263153, 263682, 269550, 269553, 269556, 269573,
                    269575, 269583, 269584, 269662, 269669, 269719, 270231, 270260,
                    276448, 276450, 276454, 449966, 463823, 569657, 600963, 600966,
                    611450, 690302, 718586, 718600, 718608, 718626, 718636, 720061,
                    720100, 754068, 754083, 768677, 803349, 803356, 807350, 808254,
                    808338, 808339, 808547, 808560, 808570, 808650, 863383]

halfed_n10_run_list = [376, 377, 378, 384, 385, 386, 798, 799, 800, 806, 807, 808, 1197, 1198, 1199, 1200, 1206, 1207]


blob_dir_dr9 = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
blob_dir_dr10 = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_ccd_blob_mask'

image_dir = '/global/project/projectdirs/cosmo/staging'
# surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-decam-dr10-v2.fits'
surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey_ccds_in_runs/survey-ccds-decam-dr10-v2-run_{}.fits'.format(run)
# output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_templates'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_templates_dev_nmad'

# ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts', 'ccdskycounts', 'plver']
ccd_columns = ['image_hdu', 'expnum', 'ccdname', 'ccdskycounts', 'ra', 'dec']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
print(len(ccd))

skyrun = Table.read('/global/cfs/cdirs/desi/users/schlafly/decals/skyrunsdr10-v2.fits')
print(len(skyrun))

mask = skyrun['run']==run
skyrun = skyrun[mask]
print(len(skyrun))

mask = skyrun['ok']==True
skyrun = skyrun[mask]
print(len(skyrun))

# # Exclude templates already created
# fn_list = glob.glob(os.path.join(output_dir, '*.fits'))
# run_list_done = [int(fn[len(os.path.join(output_dir, 'sky_templates_'))+1:-8]) for fn in fn_list]
# mask = ~np.in1d(skyrun['run'], run_list_done)
# skyrun = skyrun[mask]
# print(len(skyrun), len(run_list_done))


def get_median_stacked_image(ccdnum):

    # print(ccdnum)

    img_list = []

    ccdname = ccdnamenumdict_inv[ccdnum]

    if ccdnum<=32:
        time.sleep(1.*(ccdnum_list.index(ccdnum)))

    for index, skyrun_index in enumerate(range(len(skyrun))):

        expnum = skyrun['expnum'][skyrun_index]

        if expnum in expnum_blacklist:
            continue

        mask = ccd['expnum']==skyrun['expnum'][skyrun_index]
        ccd1 = ccd[mask].copy()

        if not ccdname in ccd1['ccdname']:
            continue

        # print(ccdnum, ccdname, index, '/', len(range(len(skyrun))))

        # Load CCD image
        img_fn = os.path.join(image_dir, skyrun['image_filename'][skyrun_index]).strip()
        ood_fn = img_fn.replace('_ooi_', '_ood_')

        if os.path.isfile(img_fn):
            img = fitsio.read(img_fn, ext=ccdname)
        else:
            print(ccdname+' '+img_fn+' does not exist!')
            continue

        if os.path.isfile(ood_fn):
            ood = fitsio.read(ood_fn, ext=ccdname)
        else:
            print(ccdname+' '+ood_fn+' does not exist!')
            continue

        if img.shape!=img_shape or img.shape!=ood.shape:
            raise ValueError

        # Get HDU index
        with fitsio.FITS(img_fn) as f:
            hdu_index = f.movnam_ext(ccdname)

        # Load blob mask
        str_loc = str.find(skyrun['image_filename'][skyrun_index].strip(), '.fits')
        img_filename_base = skyrun['image_filename'][skyrun_index].strip()[:str_loc]
        blob_path_dr9 = os.path.join(blob_dir_dr9, 'blob_mask', img_filename_base+'-blobmask.npz')
        blob_path_dr10 = os.path.join(blob_dir_dr10, 'blob_mask', img_filename_base+'-blobmask.npz')
        if os.path.isfile(blob_path_dr9) and (os.stat(blob_path_dr9).st_size!=0):
            blob_path = blob_path_dr9
        elif os.path.isfile(blob_path_dr10) and (os.stat(blob_path_dr10).st_size!=0):
            blob_path = blob_path_dr10
        else:
            # print(img_filename_base+'-blobmask.npz'+' does not exist!')
            continue
        blob_data = np.load(blob_path)

        keyname = 'hdu'+str(hdu_index).zfill(2)
        if keyname in blob_data.files:
            blob = blob_data['hdu'+str(hdu_index).zfill(2)]
        else:
            # print(blob_path+' hdu'+str(hdu_index)+' does not exist!')
            continue

        # # Find the entry in survey-ccd
        # ccd_index = np.where((ccd['expnum']==skyrun['expnum'][skyrun_index]) & (ccd['image_hdu']==hdu_index))[0][0]

        # Get the median ccdskycounts of the exposure
        # mask = (ccd['expnum']==skyrun['expnum'][skyrun_index]) & (ccd['ccdname']!='S7') & (ccd['ccdname']!='S7 ') # too slow!
        # ccdskycounts_median = np.median(ccd1['ccdskycounts'])
        # print('ccdskycounts_median = {:.4f}'.format(ccdskycounts_median))

        # Apply blob and ood mask
        img_mask = (blob==True) & (ood==0)

        # Remove median sky
        sky = np.median(img[img_mask].flatten())
        img = img - sky

        # Normalize by nmad
        img = img / nmad(img[img_mask])

        img[~img_mask] = np.nan

        # Pad to (4096, 2048) and 4x4 downsize
        img = np.pad(img, ((1, 1), (1, 1)), mode='constant', constant_values=np.nan)
        # trim_size_x = img.shape[1] % binsize
        # trim_size_y = img.shape[0] % binsize
        # img = img[:(img.shape[0]-trim_size_y), :(img.shape[1]-trim_size_x)]
        # to ignore NAN values, use np.nanmean
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore')
            img = np.nanmean(np.nanmean(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)

        img_list.append(img)

        gc.collect()

    if len(img_list)==0:
        print('WARNING: There is no available {} CCD'.format(ccdname))
        return None

    if len(img_list)<min_n_exposure:
        print('WARNING: There are only {} images available for {} CCD'.format(len(img_list), ccdname))
        # continue

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore')
        img_list = np.array(img_list)
        img_median = np.nanmedian(img_list, axis=0, overwrite_input=True)

    header = {'NIMAGE': len(img_list)}

    data = [img_median, header, ccdname]

    return data


time_start = time.time()

band = skyrun['filter'][0]

print('band: {}, run: {}'.format(band, run))

if len(skyrun)>max_n_exposure:
    np.random.seed(10000)
    idx = np.random.choice(len(skyrun), size=max_n_exposure, replace=False)
    idx.sort()
    skyrun = skyrun[idx]

output_path = os.path.join(output_dir, 'sky_raw_{}_{}.fits'.format(band, run))

if os.path.isfile(output_path):
    sys.exit(output_path+' already exists!')
    # raise ValueError(output_path+' already exists!')

with Pool(processes=n_processes) as pool:
    res = pool.map(get_median_stacked_image, ccdnum_list)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

if diagnostic_touch:
    Path('/global/u2/r/rongpu/temp/sky_template_status/'+os.path.basename(output_path)).touch()
    Path('/global/u2/r/rongpu/temp/sky_template_being_written/'+os.path.basename(output_path)).touch()

hdul_stacked = fitsio.FITS(output_path, mode='rw', clobber=True)
hdul_stacked.write(data=None)  # first HDU is empty

for tmp in res:
    hdul_stacked.write(data=tmp[0], header=tmp[1], extname=tmp[2])
hdul_stacked.close()

if diagnostic_touch:
    os.remove('/global/u2/r/rongpu/temp/sky_template_being_written/'+os.path.basename(output_path))

print('Run {} done!'.format(run), time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))

# with Pool(processes=n_processes) as pool:
#     res = pool.map(compute_raw_sky, run_list)

print('All done!!!!!!!!!!!!!!!')

