# Only process exposures with PLVER>=4.8
# Use DR9 blobs in default
# Setting n_processes to 32 might cuase out-of-memory error

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

from scipy.interpolate import RectBivariateSpline

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import search_around

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
diagnostic_plot = False
diagnostic_touch = True
overwrite = False

##########################
debug = False
##########################

image_dir = '/global/project/projectdirs/cosmo/staging/'
metrics_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/early-coadds/metrics'
coadd_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/early-coadds/coadd'

metrics_dir_dr9 = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/metrics'
coadd_dir_dr9 = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/coadd'

# surveyccd_path = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr10/survey-ccds-decam-dr10-v2.fits'
# surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey-ccds-decam-dr10-v2-trim.fits'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_ccd_blob_mask'

output_dir_dr9 = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'

# n_task = 200
# for index in range(n_task):
#     print('{} {}'.format(n_task, index))

parser = argparse.ArgumentParser()
parser.add_argument('args')
args = parser.parse_args()
n_task, task_id = args.args.split()
n_task, task_id = int(n_task), int(task_id)

surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/sky_pattern/survey_ccds_chunks/survey-ccds-decam-dr10-v2-{}.fits'.format(task_id)

# Load brick list
bricks = Table(fitsio.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr9/survey-bricks.fits.gz'))

# Load the list of CCDs left to process
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'ccdname', 'filter', 'ccd_cuts', 'ra', 'dec', 'ccdraoff', 'ccddecoff', 'plver', 'object']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))

# ###################################################
# # Only include CCDs with DEC>-25
# mask = ccd['dec']>-25
# ccd = ccd[mask]
# print('ccd', len(ccd))
# ###################################################

# mask = ccd['filter']==band
# ccd = ccd[mask]

expnum_list = np.unique(ccd['expnum'])

###################################################
if debug:
    np.random.seed(213)
    idx = np.random.choice(len(expnum_list), size=(32), replace=False)
    expnum_list = expnum_list[idx]
###################################################

print('Number of expousres in this node:', len(expnum_list))


def save_ccd_blob_mask(expnum):

    time_start = time.time()

    ccd_idx = np.where(ccd['expnum']==expnum)[0]

    str_loc = str.find(ccd['image_filename'][ccd_idx[0]].strip(), '.fits')
    img_filename_base = ccd['image_filename'][ccd_idx[0]].strip()[:str_loc]
    output_path = os.path.join(output_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
    plot_path = os.path.join(output_dir, 'plots', img_filename_base+'_blobmask_'+str(expnum))

    output_path_dr9 = os.path.join(output_dir_dr9, 'blob_mask', img_filename_base+'-blobmask.npz')
    if os.path.isfile(output_path_dr9) and (os.stat(output_path_dr9).st_size!=0):
        print('Exposure {} already exists in DR9!'.format(expnum))
        return None

    if not os.path.exists(os.path.dirname(output_path)):
        try:
            os.makedirs(os.path.dirname(output_path))
        except:
            pass

    # Check if the blobmask file already exists
    if (not overwrite) and os.path.isfile(output_path):
        print(output_path+' already exists!')
        # continue
        return None

    # Get CCD image info
    # ccd_index = len(ccd)//2
    exp_fn = os.path.join(image_dir, ccd['image_filename'][ccd_idx[0]]).strip()
    hdulist = fits.open(exp_fn)

    data = {}  # output data

    # Find bricks that cover the CCDs in this exposure
    search_radius = 1.414*((0.263*4096/2)+(0.26*3600/2))  # CCD size + brick size
    idx1, idx2, d2d, d_ra, d_dec = search_around(ccd['ra'][ccd_idx], ccd['dec'][ccd_idx], bricks['RA'], bricks['DEC'], search_radius=search_radius, verbose=False)
    mask = np.abs(d_ra)<(0.263*4096/2+0.26*3600/2)
    mask &= np.abs(d_dec)<(0.263*2048/2+0.26*3600/2)
    idx1, idx2 = idx1[mask], idx2[mask]

    all_ccds_empty = True

    for loop_index, ccd_index in enumerate(ccd_idx):

        gc.collect()

        # print(ccd_index)
        hdu_index = ccd['image_hdu'][ccd_index]
        brick_idx = idx2[idx1==loop_index]

        w_ccd = wcs.WCS(hdulist[hdu_index].header)
        naxis1 = hdulist[hdu_index].header['NAXIS1']
        naxis2 = hdulist[hdu_index].header['NAXIS2']
        # frgscale = (hdulist[hdu_index].header)['FRGSCALE']

        # here x and y are not numpy indices, but they are the x-y values in the DECam CCD schematics
        pix_x_grid, pix_y_grid = np.meshgrid(np.arange(naxis1), np.arange(naxis2))
        pix_x, pix_y = pix_x_grid.flatten(), pix_y_grid.flatten()
        # pix_ra, pix_dec = w_ccd.wcs_pix2world(pix_x, pix_y, 0)
        # Interpolate to speed up pix_ra,dec calculation:
        binsize = 100
        pix_x_spline, pix_y_spline = np.arange(-binsize, naxis1+2*binsize, binsize), np.arange(-binsize, naxis2+2*binsize, binsize)
        xx, yy = np.meshgrid(pix_x_spline, pix_y_spline)
        pix_ra_spline, pix_dec_spline = w_ccd.wcs_pix2world(xx, yy, 0)
        pix_ra_spline_wrap = ((pix_ra_spline+180.)%360.)-180.
        if not (np.all(pix_ra_spline_wrap<0) or np.all(pix_ra_spline_wrap>=0)):
            pix_ra_spline = pix_ra_spline_wrap
        pix_ra_spline = pix_ra_spline + ccd['ccdraoff'][ccd_index] / 3600.
        pix_dec_spline = pix_dec_spline + ccd['ccddecoff'][ccd_index] / 3600.

        # # Find bricks that cover the CCD --- OLD VERSION THAT DOES NOT ALWAYS WORK!!!
        # ccd_corners = [[pix_ra.min(), pix_dec.min()],
        #                [pix_ra.min(), pix_dec.max()],
        #                [pix_ra.max(), pix_dec.min()],
        #                [pix_ra.max(), pix_dec.max()]]
        # ccd_corners = np.array(ccd_corners)
        # mask = np.zeros(len(bricks), dtype=bool)
        # for index in range(len(ccd_corners)):
        #     mask |= (bricks['RA1']<ccd_corners[index, 0]) & (bricks['RA2']>ccd_corners[index, 0]) \
        #      & (bricks['DEC1']<ccd_corners[index, 1]) & (bricks['DEC2']>ccd_corners[index, 1])
        # print(np.sum(mask))
        # brick_idx = np.where(mask)[0]

        # Further reduce the number of non-overlapping bricks
        # Only DEC limits; RA limits would be more complicated
        mask = (pix_dec_spline.min()<=bricks['DEC2'][brick_idx]) & (pix_dec_spline.max()>=bricks['DEC1'][brick_idx])
        brick_idx = brick_idx[mask]
        # print('expnum={}, loop={}, {} bricks'.format(expnum, loop_index, len(brick_idx)))

        # Diagnostic plot to check the brick coverage
        if diagnostic_plot:
            if not os.path.exists(os.path.dirname(plot_path)):
                try:
                    os.makedirs(os.path.dirname(plot_path))
                except:
                    pass
            points = w_ccd.wcs_pix2world([[0, 0], [naxis1, naxis2], [0, naxis2], [naxis1, 0]], 0)
            plt.figure(figsize=(6, 6))
            plt.plot(pix_ra_spline, pix_dec_spline, '.', ms=0.5)
            for index in brick_idx:
                ra_plot = [bricks['RA1'][index], bricks['RA1'][index], bricks['RA2'][index], bricks['RA2'][index], bricks['RA1'][index]]
                dec_plot = [bricks['DEC1'][index], bricks['DEC2'][index], bricks['DEC2'][index], bricks['DEC1'][index], bricks['DEC1'][index]]
                plt.plot(ra_plot, dec_plot, alpha=0.5)
            plt.plot(points[:, 0], points[:, 1], '.', ms=5, color='r')
            plt.xlabel('RA')
            plt.ylabel('DEC')
            # plt.tight_layout()
            plt.savefig(plot_path+'_hdu{}_1.png'.format(str(hdu_index).zfill(2)))
            plt.close()

        # Obtain the blob mask
        img_mask = np.zeros([naxis2, naxis1], dtype=bool)
        for brick_index in brick_idx:
            brickname = bricks['BRICKNAME'][brick_index]
            # ra1, ra2, dec1, dec2 = bricks['RA1'][index], bricks['RA2'][index], bricks['DEC1'][index], bricks['DEC2'][index]
            # there are bricks are missing
            blobs_path = os.path.join(metrics_dir, '{}/blobmask-{}.fits.gz'.format(brickname[:3], brickname))
            maskbits_path = os.path.join(coadd_dir, '{}/{}/legacysurvey-{}-maskbits-light.fits.fz'.format(brickname[:3], brickname, brickname))
            blobs_path_dr9 = os.path.join(metrics_dir_dr9, '{}/blobs-{}.fits.gz'.format(brickname[:3], brickname))
            maskbits_path_dr9 = os.path.join(coadd_dir_dr9, '{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(brickname[:3], brickname, brickname))

            if os.path.isfile(blobs_path_dr9):  # Use DR9 blobs by default
                blobs = fits.getdata(blobs_path_dr9)
                maskbits = fits.getdata(maskbits_path_dr9)
                blob_hdu = fits.open(blobs_path_dr9)
                mask_bad = (blobs!=-1)  # -1 means no source detected
                mask_bad |= (maskbits&2**0>0)  # brick_primary
                mask_bad |= (maskbits&2**1>0)  # BRIGHT star
                mask_bad |= (maskbits&2**5>0) | (maskbits&2**6>0) | (maskbits&2**7>0)  # allmask
                mask_bad |= (maskbits&2**11>0) | (maskbits&2**12>0) | (maskbits&2**13>0)  # MEDIUM, GALAXY, CLUSTER
            elif (not 'DES supernova' in ccd['object'][ccd_index]) and os.path.isfile(blobs_path):
                blobs = fits.getdata(blobs_path)
                maskbits = fits.getdata(maskbits_path)
                blob_hdu = fits.open(blobs_path)
                mask_bad = (blobs==1)
                mask_bad |= (maskbits&2**0>0)  # brick_primary
                mask_bad |= (maskbits&2**1>0)  # BRIGHT star
                mask_bad |= (maskbits&2**11>0) | (maskbits&2**12>0) | (maskbits&2**13>0)  # MEDIUM, GALAXY, CLUSTER
            else:
                continue
            mask_good = ~mask_bad

            w_blob = wcs.WCS(blob_hdu[0].header)
            coadd_x_spline, coadd_y_spline = w_blob.wcs_world2pix(pix_ra_spline, pix_dec_spline, 0)
            interp_x = RectBivariateSpline(pix_y_spline, pix_x_spline, coadd_x_spline)
            interp_y = RectBivariateSpline(pix_y_spline, pix_x_spline, coadd_y_spline)
            coadd_x = interp_x(np.arange(naxis2), np.arange(naxis1)).flatten()
            coadd_y = interp_y(np.arange(naxis2), np.arange(naxis1)).flatten()
            coadd_x, coadd_y = np.round(coadd_x).astype(int), np.round(coadd_y).astype(int)
            mask = (coadd_x>=0) & (coadd_x<blobs.shape[0]) & (coadd_y>=0) & (coadd_y<blobs.shape[0])
            img_mask[pix_y[mask], pix_x[mask]] |= mask_good[coadd_y[mask], coadd_x[mask]]

        if np.sum(img_mask)!=0:
            all_ccds_empty = False
            data['hdu'+str(hdu_index).zfill(2)] = img_mask

        # Diagnostic plot to check the blobmask image
        if diagnostic_plot:
            fig = plt.figure(frameon=False, figsize=(8, 4))
            ax = fig.add_axes([0, 0, 1, 1])
            ax.axis('off')
            ax.imshow(img_mask.T, cmap='gray', origin='upper')
            # plt.tight_layout()
            plt.savefig(plot_path+'_hdu{}_2.png'.format(str(hdu_index).zfill(2)))
            plt.close()

    #################################################################################
    # Checking again because another process might running on the same exposures
    if (not overwrite) and os.path.isfile(output_path):
        print(output_path+' already exists!')
        # continue
        return None
    #################################################################################

    # Save blob mask image
    if not all_ccds_empty:
        if diagnostic_touch:
            # Path('/global/cfs/cdirs/desi/users/rongpu/tmp/blobmask_status/'+os.path.basename(output_path)).touch()
            Path('/global/cfs/cdirs/desi/users/rongpu/tmp/blobmask_being_written/'+os.path.basename(output_path)).touch()
        np.savez_compressed(output_path, **data)
        if diagnostic_touch:
            os.remove('/global/cfs/cdirs/desi/users/rongpu/tmp/blobmask_being_written/'+os.path.basename(output_path))
        print('Exposure {} done!'.format(expnum), time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
    else:
        print('Exposure {} is not in the footprint!!!'.format(expnum))
        Path(output_path).touch()

    return None


# def pool_wrapper(expnum_list):
#     if len(expnum_list)==0:
#         return None
#     for expnum in expnum_list:
#         save_ccd_blob_mask(expnum)
#     return None

##############################################################################################################################

def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(save_ccd_blob_mask, expnum_list, chunksize=1)

    print('All done!!!!!!!!!!!!!!!')


if __name__=="__main__":
    main()
