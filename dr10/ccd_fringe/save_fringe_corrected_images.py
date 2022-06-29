from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from multiprocessing import Pool
from pathlib import Path

params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

time_start_all = time.time()

n_processes = 64

keep_status = True

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

fringe_dr9_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/fringe/DECam_CP-Fringe-Normed'
fringe_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/fringe_templates'

image_dir = '/global/cfs/cdirs/cosmo/staging/'
surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1.fits'

blob_dir_dr9 = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
blob_dir_dr10 = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_ccd_blob_mask'

frgscale_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/frgscale'

image_output_dir = '/pscratch/sd/r/rongpu/dr10dev/fringe_corrected_images'

status_dir = '/global/homes/r/rongpu/temp/status'


##############################################################################################################################

# Load CCD list
# ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts', 'mjd_obs']
# ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
mask = ccd['filter']=='z'  # include only z-band images
ccd = ccd[mask]
print(len(ccd), len(np.unique(ccd['expnum'])))
ccd['ccdnum'] = [ccdnamenumdict[ccd['ccdname'][ii].strip()] for ii in range(len(ccd))]

expnum_list = np.unique(ccd['expnum'])

# # Remove exposures that are already done
# expnum_list_done = np.zeros(len(expnum_list), dtype=bool)
# for index, expnum in enumerate(expnum_list):
#     # Find an arbitrary CCD in the exposure to get the image filename
#     ccd_index = np.where((ccd['expnum']==expnum))[0][0]
#     img_fn_write = os.path.join(image_output_dir, ccd['image_filename'][ccd_index].strip())
#     if os.path.isfile(img_fn_write):
#         expnum_list_done[index] = True
# print('Done     Not-done    Done/Not-done')
# print(np.sum(expnum_list_done), np.sum(~expnum_list_done), np.sum(expnum_list_done)/len(expnum_list_done))
# expnum_list = expnum_list[~expnum_list_done]
# print('Expsoures left to process: ', len(expnum_list))

# shuffle
np.random.seed(12345)
# DO NOT USE NP.RANDOM.SHUFFLE
expnum_list = np.random.choice(expnum_list, size=len(expnum_list), replace=False)

# Load DR9 fringe templates
fringe_templates_dr9 = {}
for ccdnum in range(1, 63):
    ccdname = ccdnamenumdict_inv[ccdnum]
    fringe_template_path = os.path.join(fringe_dr9_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
    if os.path.isfile(fringe_template_path):
        fringe_tmp = fitsio.read(fringe_template_path)
        fringe_tmp = fringe_tmp[1:4095, 1:2047]  # remove the edge pixels
        fringe_templates_dr9[ccdnum] = fringe_tmp.copy()

# Load initial fringe templates
fringe_templates = {}
for ccdnum in range(1, 63):
    ccdname = ccdnamenumdict_inv[ccdnum]
    fringe_template_path = os.path.join(fringe_dir, 'DECam_z_frg_{}_CCD{}.fits'.format(ccdname, str(ccdnum).zfill(2)))
    if os.path.isfile(fringe_template_path):
        fringe_tmp = fitsio.read(fringe_template_path)
        fringe_tmp = fringe_tmp[1:4095, 1:2047]  # remove the edge pixels
        fringe_templates[ccdnum] = fringe_tmp.copy()
print('fringe_templates', len(fringe_templates))

fringe_table = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/frgscale/survey-ccds-dr10-deep-fields-v1-frgscales.fits'))
fringe_table_all = fringe_table.copy()

fringe_stats = Table.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/frgscale/survey-ccds-dr10-deep-fields-v1-frgscales-stats.txt', format='ascii.commented_header')

##############################################################################################################################


def save_image(expnum):

    time_start = time.time()

    print('expnum:', expnum)

    # Find an arbitrary CCD in the exposure to get the image filename
    ccd_index = np.where((ccd['expnum']==expnum))[0][0]

    # frgscale_path = os.path.join(frgscale_dir, ccd['image_filename'][ccd_index].strip().replace('.fits.fz', '.txt'))
    img_fn_write = os.path.join(image_output_dir, ccd['image_filename'][ccd_index].strip())
    img_fn = os.path.join(image_dir, ccd['image_filename'][ccd_index].strip())
    print(img_fn_write)

    if os.path.isfile(img_fn_write):
        print(img_fn_write+' already exists!')
        return None

    if keep_status:
        status_fn = os.path.join(status_dir, str(expnum))
        Path(status_fn).touch()

    if not os.path.exists(os.path.dirname(img_fn_write)):
        try:
            os.makedirs(os.path.dirname(img_fn_write))
        except:
            pass

    # if not os.path.isfile(frgscale_path):
    #     print('Frgscale does not exist', frgscale_path)
    #     return None
    # elif (os.stat(frgscale_path).st_size==0):
    #     print('Frgscale is empty', frgscale_path)
    #     return None
    # # fringe_table = Table.read(frgscale_path, format='ascii.commented_header')

    mask = fringe_table_all['expnum']==expnum
    fringe_table = fringe_table_all[mask]
    mask = (fringe_table['ccdname']!='S7') & (fringe_table['ccdname']!=['S30'])
    n_ccd = np.sum(mask)

    if n_ccd<10:
        print('Only {} CCDs available'.format(n_ccd))

    hdul_r = fitsio.FITS(img_fn, mode='r')
    hdul_w = fitsio.FITS(img_fn_write, mode='rw', clobber=True)

    for hdu_index in range(len(hdul_r)):
        if hdu_index==0:
            hdr = hdul_r[hdu_index].read_header()
            hdul_w.write(data=None, header=hdr)
        else:
            hdr = hdul_r[hdu_index].read_header()
            img = hdul_r[hdu_index].read()
            ccdname = hdul_r[hdu_index].get_extname()
            ccdnum = ccdnamenumdict[ccdname.strip()]

            # Back out the exisiting fringe correction
            if 'FRGSCNEW' in hdr:
                frgscale_dr9 = hdr['FRGSCNEW']
                img += fringe_templates_dr9[ccdnum] * frgscale_dr9

            mask = (fringe_table['ccdname']==ccdname)
            if np.sum(mask)==1:
                frgscale = fringe_table['frgscale_apply'][mask][0]
            else:
                fringe_stats_index = np.where(fringe_stats['ccdname']==ccdname)[0][0]
                frgscale = fringe_table['frgscale_median'][0] * fringe_stats['median_ratio'][fringe_stats_index]

            img -= frgscale * fringe_templates[ccdnum]
            new_key = {'name': 'FRGSCLV2', 'value': frgscale, 'comment': 'New (v2) fringe-correction scale'}
            hdr.add_record(new_key)

            hdul_w.write(data=img, header=hdr, extname=ccdname, compress='rice')

    hdul_r.close()
    hdul_w.close()

    gc.collect()

    if keep_status:
        os.remove(status_fn)

    print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))

    return None


with Pool(processes=n_processes) as pool:
    res = pool.map(save_image, expnum_list)

print('All!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start_all)))


