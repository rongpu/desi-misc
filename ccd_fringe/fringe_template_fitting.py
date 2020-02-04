from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits
from astropy import wcs

from multiprocessing import Pool
import argparse
from pathlib import Path

from scipy.interpolate import RectBivariateSpline
from scipy.ndimage.filters import gaussian_filter
from scipy import stats

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'}
plt.rcParams.update(params)

n_node = 2 # Haswell
n_processess = 31

parser = argparse.ArgumentParser()
parser.add_argument('task_id')
args = parser.parse_args()
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

fringe_old_dir = '/global/homes/d/djschleg/cosmo/staging/decam/DECam_CP-Fringe'
fringe_new_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/DECam_CP-Fringe'
image_dir = '/global/project/projectdirs/cosmo/staging/'
# surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'
surveyccd_path = '/global/homes/r/rongpu/mydesi/dr9/fringe/misc/survey-ccds-decam-dr9-z-band-only-trim.fits'
blob_dir = '/global/cscratch1/sd/rongpu/fringe/decam_ccd_blob_mask'

# image_output_dir = '/global/cscratch1/sd/rongpu/fringe/fringe_corrected_image/'
frgscale_output_dir = '/global/cscratch1/sd/rongpu/fringe/frgscale/'
# image_output_dir = '/global/cscratch1/sd/rongpu/fringe/tmp_img/'
# frgscale_output_dir = '/global/cscratch1/sd/rongpu/fringe/tmp_frgscale/'

# expnum_list = np.array([243575, 257579, 608237, 613691, 615457, 617149, 625647, 625963, 625984,
#  626013, 626022, 626369, 626404, 630676, 630764, 630767, 630920, 630975,
#  631032, 631079, 635121, 637641, 640631, 640713, 641990, 642213, 648115,
#  648436, 648438, 648446, 648449, 648503, 649869, 649876, 649925, 660071,
#  662042, 663660, 666027, 666048, 675281, 675757, 676897, 677354, 677384,
#  678507, 683976, 685495, 685846, 685852, 690311, 690318, 690446, 690452,
#  690465, 690480, 690495, 690809, 690810, 691185, 692345, 693115, 693133,
#  695083, 695086, 695092, 695470, 695823, 696082, 698757, 698762, 699867,
#  700267, 702094, 702347, 705161, 714405, 715780, 716059, 716111, 717412,
#  718693, 719083, 719086, 719104, 719810, 720123, 720132, 720140, 721256,
#  721400, 723171, 723176, 724737, 725084, 725092, 725112, 730384, 731280,
#  731283, 731401, 745654, 746829, 747025, 759175, 763430, 764325, 764460,
#  766535, 766609, 766944, 767015, 767038, 767350, 767357, 767359, 767361,
#  767551, 767553, 767576, 767946, 767975, 768022, 768430, 768461, 768471,
#  768829, 769445, 770830, 773167, 773187, 773190, 773232, 773628, 773643,
#  774102, 774145, 774226, 774566, 774576, 774644, 775024, 775050, 775054,
#  775090, 775921, 775959, 776279, 776305, 776310, 779802, 779935, 779962,
#  780012, 780257, 780282, 781167, 781516, 783015, 783532, 783557, 783558,
#  783576, 783821, 783823, 789677, 790034, 790358, 792830, 793321, 793324,
#  793349, 802014, 807112, 807394, 807403, 829963, 830087, 830350, 830662])
# print(len(expnum_list))

##############################################################################################################################

# Trim CCD edges
img_trim_size_x = 23 + 200
img_trim_size_y = 47 + 200

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts']
ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
# ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
mask = ccd['ccd_cuts']==0
mask &= ccd['filter']=='z' # include only z-band images
ccd = ccd[mask]
print(len(ccd))
ccd['ccdnum'] = [ccdnamenumdict[ccd['ccdname'][ii].strip()] for ii in range(len(ccd))]

expnum_list = np.unique(ccd['expnum'])
# shuffle
np.random.seed(123)
# DO NOT USE NP.RANDOM.SHUFFLE
expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)

# split among the Cori nodes
expnum_list_split = np.array_split(expnum_list, n_node)
expnum_list = expnum_list_split[task_id]


# Load old fringe images
fringe_old_dict = {}
for ccdnum in range(1, 63):
    # skip N30 and S7
    if ccdnum==61 or ccdnum==31:
        continue
    fringe_old_path = os.path.join(fringe_old_dir, 'DES17B_20180103_908c062-z-{}_frg.fits'.format(str(ccdnum).zfill(2)))
    fringe_old = fits.getdata(fringe_old_path)
    # remove the edge pixels
    fringe_old = fringe_old[1:4095, 1:2047]
    fringe_old_dict[ccdnum] = fringe_old.copy()

# Load new fringe images
fringe_dict_original = {}
fringe_dict_spline_subtracted = {}
for ccdnum in range(1, 63):
    # skip N30 and S7
    if ccdnum==61 or ccdnum==31:
        continue
    fringe_path = glob.glob(os.path.join(fringe_new_dir, '*CCD{}.fits'.format(str(ccdnum).zfill(2))))[0]
    fringe = fits.getdata(fringe_path)
    fringe = fringe[1:4095, 1:2047]
    fringe_dict_original[ccdnum] = fringe.copy()
    
    # Trim CCD edges
    fringe = fringe[img_trim_size_y:(fringe.shape[0]-img_trim_size_y), img_trim_size_x:(fringe.shape[1]-img_trim_size_x)]

    # Spline sky: fringe
    # downsize the image to speed up computation
    binsize = 400
    trim_size_x = 0
    trim_size_y = 0
    fringe_spline_data = np.nanmedian(np.nanmedian(fringe.reshape((fringe.shape[0]//binsize, binsize, fringe.shape[1]//binsize,-1)), axis=3), axis=1)
    x_sky_grid = trim_size_x + binsize/2+binsize*np.arange(fringe_spline_data.shape[1])
    y_sky_grid = trim_size_y + binsize/2+binsize*np.arange(fringe_spline_data.shape[0])
    # if np.sum(~np.isfinite(fringe_spline_data))!=0:
        # print(np.sum(~np.isfinite(fringe_spline_data)), 'NAN values!')
    spline = RectBivariateSpline(y_sky_grid, x_sky_grid, fringe_spline_data)
    fringe_spline = spline(np.arange(fringe.shape[0]), np.arange(fringe.shape[1]))
    # Subtract spline
    fringe -= fringe_spline
    fringe_dict_spline_subtracted[ccdnum] = fringe.copy()

##############################################################################################################################

# expnum = 718884
# expnum = 624440
# expnum = expnum_list[0]

def compute_fringe_scale(expnum):

    print('expnum:', expnum)

    # Find an arbitrary CCD in the exposure to get the image filename
    ccd_index = np.where((ccd['expnum']==expnum))[0][0]

    frgscale_output_path = os.path.join(frgscale_output_dir, ccd['image_filename'][ccd_index].strip().replace('.fits.fz', '.txt'))
    if os.path.isfile(frgscale_output_path):
        print(frgscale_output_path+' already exists!')
        return None

    if not os.path.exists(os.path.dirname(frgscale_output_path)):
        try:
            os.makedirs(os.path.dirname(frgscale_output_path))
        except:
            pass

    # Load blob mask
    blob_path = os.path.join(blob_dir, 'blob_mask', ccd['image_filename'][ccd_index].strip().replace('.fits.fz', '-blobmask.npz'))
    
    try:
        blob_data = np.load(blob_path)
    except:
        print('Exposure {} is not in the DR8 footprint!!!'.format(expnum))
        Path(frgscale_output_path).touch()
        return None

    img_fn = os.path.join(image_dir, ccd['image_filename'][ccd_index].strip())
    hdulist = fits.open(img_fn)

    fringe_params = []

    for hdu_index in range(1, len(hdulist)):
    # for hdu_index in [3]:

        hdulist = fits.open(img_fn)

        # CP mask image
        ood_fn = img_fn.replace('_ooi_', '_ood_')
        ood_hdulist = fits.open(ood_fn)
        
        ccdname = hdulist[hdu_index].header['EXTNAME'].strip()
        ccdnum = ccdnamenumdict[ccdname]

        # Some images do not have FRGSCALE in the header
        # (they were not fringed corrected in CP, so we skip them here as well)
        try:
            frgscale = (hdulist[hdu_index].header)['FRGSCALE']
        except:
            if ccdname!='S7':
                print('Error: no frgscale for {}'.format(ccdname))
            # else:
            #     print('no frgscale for {}'.format(ccdname))
            continue

        # print('HDU:', hdu_index)
            
        fringe = fringe_dict_spline_subtracted[ccdnum]
        fringe_original = fringe_dict_original[ccdnum]

        # Load CCD image
        img = hdulist[hdu_index].data.copy()

        # Back out the exisiting fringe correction
        fringe_old = fringe_old_dict[ccdnum]
        img += fringe_old*frgscale
        img_original = img.copy()
        
        # CP mask
        ood = ood_hdulist[hdu_index].data

        # Blob mask
        try:
            blob = blob_data['hdu'+str(hdu_index).zfill(2)]
        except:  # happens when the CCD is outside the DR8 blobmask coverage
            continue

        # Apply masks
        # blob: no source is True
        img_mask = (blob==True) & (ood==0)
        if np.sum(img_mask)==0:
            continue

        img[~img_mask] = np.nan

        # Remove median sky
        median_sky = np.median(img[np.isfinite(img)])
        img = img - median_sky

        # # Normalize by frgscale
        # img = img/frgscale

        # Trim CCD edges
        img = img[img_trim_size_y:(img.shape[0]-img_trim_size_y), img_trim_size_x:(img.shape[1]-img_trim_size_x)]
        # print(img.shape)

        # 3-sigma clipping
        sky_nmad = nmad(img[np.isfinite(img)]) # sky level
        # print('sky nmad:', sky_nmad)
        mask = (img<-3*sky_nmad) | (img>3*sky_nmad)
        img[mask] = 0
        
        ##################################################################################################################

        # Spline sky: image
        # downsize the image to speed up computation
        binsize = 400
        trim_size_x = 0
        trim_size_y = 0
        img_spline_data = np.nanmedian(np.nanmedian(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)
        x_sky_grid = trim_size_x + binsize/2+binsize*np.arange(img_spline_data.shape[1])
        y_sky_grid = trim_size_y + binsize/2+binsize*np.arange(img_spline_data.shape[0])
        if np.sum(~np.isfinite(img_spline_data))!=0:
            print(ccdnum, np.sum(~np.isfinite(img_spline_data)), 'NAN values!')
        mask = ~np.isfinite(img_spline_data)
        img_spline_data[mask] = 0

        spline = RectBivariateSpline(y_sky_grid, x_sky_grid, img_spline_data)
        img_spline = spline(np.arange(img.shape[0]), np.arange(img.shape[1]))

        # Subtract spline
        img -= img_spline

        # Linear regression
        img_mask = np.isfinite(img)
        if np.sum(img_mask)==0:
            continue
        img_flat = img[img_mask].flatten()
        fringe_flat = fringe[img_mask].flatten()

        slope, intercept, r_value, p_value, std_err = stats.linregress(fringe_flat, img_flat)
        # print(slope, intercept, r_value, p_value, std_err)
        
        fringe_params.append([expnum, ccdnum, frgscale, median_sky, sky_nmad, slope, intercept, r_value, p_value, std_err])

        # hdulist[hdu_index].data = img_original - slope * fringe_original
        
        # img1 = img.copy()
        # img1[~np.isfinite(img1)] = 0
        # img1 = gaussian_filter(img1, 4, mode='reflect', truncate=3)
        # vrange = 0.5*sky_nmad
        # plt.figure(figsize=(10, 5))
        # plt.imshow((img1).T, cmap='seismic', vmin=-vrange, vmax=vrange)
        # plt.tight_layout()
        # plt.show()

        # img1 = img - slope*fringe
        # img1[~np.isfinite(img1)] = 0
        # img1 = gaussian_filter(img1, 4, mode='reflect', truncate=3)
        # vrange = 0.5*sky_nmad
        # plt.figure(figsize=(10, 5))
        # plt.imshow((img1).T, cmap='seismic', vmin=-vrange, vmax=vrange)
        # plt.tight_layout()
        # plt.show()

        hdulist.close()
        ood_hdulist.close()

    if len(fringe_params)==0:
        Path(frgscale_output_path).touch()
    else:        
        fringe_table = Table(np.array(fringe_params), names=['expnum', 'ccdnum', 'frgscale_old', 'median_sky', 'sky_nmad', 'slope', 'intercept', 'r_value', 'p_value', 'std_err'])
        fringe_table['expnum'] = np.array(fringe_table['expnum'], dtype=int)
        fringe_table['ccdnum'] = np.array(fringe_table['ccdnum'], dtype=int)
        fringe_table.write(frgscale_output_path, format='ascii.commented_header', overwrite=True)


def main():

    with Pool(processes=n_processess) as pool:
        res = pool.map(compute_fringe_scale, expnum_list)

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

