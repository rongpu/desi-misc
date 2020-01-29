# Based on get_blob_mask_for_ccd

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
sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import search_around

from multiprocessing import Pool

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'}
plt.rcParams.update(params)

n_processess = 60
diagnostic_plot = True

image_dir = '/global/project/projectdirs/cosmo/staging/'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9-cut.fits.gz'
output_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
# output_dir = '/global/homes/r/rongpu/data/decam_ccd_blob_mask'

##############################################################################################################################

# Load brick list
bricks = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/survey-bricks.fits.gz')

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'filter', 'mjd_obs', 'ra', 'dec', 'ccdraoff', 'ccddecoff', 'ccd_cuts']
ccd = fitsio.read(surveyccd_path, columns=ccd_columns)
# ccd = fitsio.read(surveyccd_path)
ccd = Table(ccd)
mask = ccd['ccd_cuts']==0
mask &= ccd['filter']=='z' # include only z-band images
ccd = ccd[mask]
print(len(ccd), 'CCD')

# # Find CCDs from a particular date
# mask = np.char.find(np.array(ccd['image_filename'], dtype='str'), 'CP20170304')!=-1
# ccd = ccd[mask]
# print(len(ccd), 'CCDs')
# # Remove the dawn frames
# mjd_sort = np.sort(ccd['mjd_obs'])
# mdj_sep = (mjd_sort[1:] - mjd_sort[:-1])
# mask = ccd['mjd_obs']<=mjd_sort[np.argmax(mdj_sep)]
# ccd = ccd[mask]
# print(len(ccd))

# Find CCDs around some MJD
mask = (ccd['mjd_obs']>(57815-4)) & (ccd['mjd_obs']<(57815+4)) # DECaLS observing run starting Feb 28, 2017
mask |= ((ccd['mjd_obs']>(58359-2)) & (ccd['mjd_obs']<(58359+27))) # Starting Aug 28, 2018
mask |= ((ccd['mjd_obs']>(58423-2)) & (ccd['mjd_obs']<(58423+30))) # Two runs starting Oct 28, 2018
mask |= ((ccd['mjd_obs']>(57893-2)) & (ccd['mjd_obs']<(57893+30))) # Two runs starting May 18, 2017
ccd = ccd[mask]
print(len(ccd), 'CCDs')

expnum_list = np.unique(ccd['expnum'])
print('Total nubmer of exposures:', len(expnum_list))

# Find the CCDs whose blobmask files do not yet exist
ccd_mask = np.zeros(len(ccd), dtype=bool) # True if exist
for ccd_index in range(len(ccd)):
    str_loc = str.find(ccd['image_filename'][ccd_index].strip(), '.fits')
    img_filename_base = ccd['image_filename'][ccd_index].strip()[:str_loc]
    blob_path = os.path.join(output_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
    if os.path.isfile(blob_path):
        ccd_mask[ccd_index] = True
print(np.sum(ccd_mask)/len(ccd_mask))
ccd = ccd[~ccd_mask]
print(len(ccd))
expnum_list = np.unique(ccd['expnum'])
expnum_list.sort()
print('Nubmer of exposures left to process:', len(expnum_list))
print()

##############################################################################################################################

def save_ccd_blob_mask(expnum):

    ccd_idx = np.where(ccd['expnum']==expnum)[0]

    str_loc = str.find(ccd['image_filename'][ccd_idx[0]].strip(), '.fits')
    img_filename_base = ccd['image_filename'][ccd_idx[0]].strip()[:str_loc]
    output_path = os.path.join(output_dir, 'blob_mask', img_filename_base+'-blobmask.npz')
    plot_path = os.path.join(output_dir, 'plots', img_filename_base+'-blobmask')
    if not os.path.exists(os.path.dirname(plot_path)):
        os.makedirs(os.path.dirname(plot_path))
    if not os.path.exists(os.path.dirname(output_path)):
        os.makedirs(os.path.dirname(output_path))

    # Check if the blobmask file already exists
    if os.path.isfile(output_path):
        print(output_path+' already exists!')
        # continue
        return None

    # Get CCD image info
    # ccd_index = len(ccd)//2
    exp_fn = os.path.join(image_dir, ccd['image_filename'][ccd_idx[0]]).strip()
    hdulist = fits.open(exp_fn)

    data = {} # output data

    # Find bricks that cover the CCDs in this exposure
    search_radius = 1.414*((0.263*4096/2)+(0.26*3600/2)) # CCD size + brick size
    idx1, idx2, d2d, d_ra, d_dec = search_around(ccd['ra'][ccd_idx], ccd['dec'][ccd_idx], bricks['RA'], bricks['DEC'], search_radius=search_radius, verbose=False)
    mask = np.abs(d_ra)<(0.263*4096/2+0.26*3600/2)
    mask &= np.abs(d_dec)<(0.263*2048/2+0.26*3600/2)
    idx1, idx2 = idx1[mask], idx2[mask]

    for loop_index, ccd_index in enumerate(ccd_idx):

        # print(ccd_index)
        hdu_index = ccd['image_hdu'][ccd_index]
        brick_idx = idx2[idx1==loop_index]

        w = wcs.WCS(hdulist[hdu_index].header)
        naxis1 = hdulist[hdu_index].header['NAXIS1']
        naxis2 = hdulist[hdu_index].header['NAXIS2']
        # frgscale = (hdulist[hdu_index].header)['FRGSCALE']

        # here x and y are not numpy indices, but they are the x-y values in the DECam CCD schematics
        pix_x_grid, pix_y_grid = np.meshgrid(np.arange(naxis1), np.arange(naxis2))
        pix_x, pix_y = pix_x_grid.flatten(), pix_y_grid.flatten()
        pix_ra, pix_dec = w.wcs_pix2world(pix_x, pix_y, 0)
        pix_ra = pix_ra + ccd['ccdraoff'][ccd_index] / 3600.
        pix_dec = pix_dec + ccd['ccddecoff'][ccd_index] / 3600.

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
        mask = (pix_dec.min()<=bricks['DEC2'][brick_idx]) & (pix_dec.max()>=bricks['DEC1'][brick_idx])
        brick_idx = brick_idx[mask]
        print('expnum={}, loop={}, {} bricks'.format(expnum, loop_index, len(brick_idx)))

        # Diagnostic plot to check the brick coverage
        if diagnostic_plot:
            points = w.wcs_pix2world([[0,0], [naxis1, naxis2], [0, naxis2], [naxis1, 0]], 0)
            plt.figure(figsize=(6, 6))
            plt.plot(pix_ra[::500], pix_dec[::500], '.', ms=0.5)
            for index in brick_idx:
                ra_plot = [bricks['RA1'][index], bricks['RA1'][index], bricks['RA2'][index], bricks['RA2'][index], bricks['RA1'][index]]
                dec_plot = [bricks['DEC1'][index], bricks['DEC2'][index], bricks['DEC2'][index], bricks['DEC1'][index], bricks['DEC1'][index]]
                plt.plot(ra_plot, dec_plot, alpha=0.5)
            plt.plot(points[:, 0], points[:, 1], '.', ms=5, color='r')
            plt.xlabel('RA')
            plt.ylabel('DEC')
            plt.savefig(plot_path+'_hdu{}_1.png'.format(str(hdu_index).zfill(2)))
            plt.close()

        # Obtain the blob mask
        img_mask = np.zeros([naxis2, naxis1], dtype=bool)
        for brick_index in brick_idx:
            brickname = bricks['BRICKNAME'][brick_index]
            # there are bricks are missing
            try:
                blobs = fits.getdata('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/south/metrics/{}/blobs-{}.fits.gz'.format(brickname[:3], brickname))
                maskbits = fits.getdata('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/south/coadd/{}/{}/legacysurvey-{}-maskbits.fits.fz'.format(brickname[:3], brickname, brickname))
            except:
                print('Warning: brick {} does not exist'.format(brickname))
                continue
            mask_bad = (blobs!=-1) # -1 means no source detected
            mask_bad |= (maskbits&2**0>0) # brick_primary
            mask_bad |= (maskbits&2**1>0) # BRIGHT star
            mask_bad |= (maskbits&2**2>0) | (maskbits&2**3>0) | (maskbits&2**4>0)  # saturation
            mask_bad |= (maskbits&2**5>0) | (maskbits&2**6>0) | (maskbits&2**7>0)  # allmask
            mask_bad |= (maskbits&2**11>0) | (maskbits&2**12>0) | (maskbits&2**13>0) # MEDIUM, GALAXY, CLUSTER
            mask_good = ~mask_bad
            blob_hdu = fits.open('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/south/metrics/{}/blobs-{}.fits.gz'.format(brickname[:3], brickname))
            w_blob = wcs.WCS(blob_hdu[0].header)
            coadd_x, coadd_y = w_blob.wcs_world2pix(pix_ra, pix_dec, 0)
            coadd_x, coadd_y = np.round(coadd_x).astype(int), np.round(coadd_y).astype(int)
            mask = (coadd_x>=0) & (coadd_x<blobs.shape[0]) & (coadd_y>=0) & (coadd_y<blobs.shape[0])
            img_mask[pix_y[mask], pix_x[mask]] |= mask_good[coadd_y[mask], coadd_x[mask]]

        # Diagnostic plot to check the blobmask image
        if diagnostic_plot:
            fig = plt.figure(frameon=False, figsize=(8, 4))
            ax = fig.add_axes([0, 0, 1, 1])
            ax.axis('off')
            ax.imshow(img_mask.T, cmap='gray', origin='upper')
            plt.savefig(plot_path+'_hdu{}_2.png'.format(str(hdu_index).zfill(2)))
            plt.close()

        data['hdu'+str(hdu_index).zfill(2)] = img_mask

    # Save blob mask image
    np.savez_compressed(output_path, **data)

    return None

def pool_wrapper(expnum_list):
    if len(expnum_list)==0:
        return None
    for expnum in expnum_list:
        save_ccd_blob_mask(expnum)
    return None

##############################################################################################################################

def main():
    
    # # only process a single exposure
    # save_ccd_blob_mask(expnum_list[0])

    # parralism via multiprocessing
    # expnum_list_split = np.array_split(expnum_list, n_processess)
    idx = np.arange(len(expnum_list))
    expnum_list_split = []
    for index in range(n_processess):
        indices = idx[::n_processess]+index
        indices = indices[indices<len(expnum_list)]
        expnum_list_split.append(expnum_list[indices])
    print()

    with Pool(processes=n_processess) as pool:
        res = pool.map(pool_wrapper, expnum_list_split)

if __name__=="__main__":
    main()
