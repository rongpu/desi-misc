# Find overlapping bricks of 90prime CCDs

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

n_processes = 32

pixel_size = 0.454

image_dir = '/global/project/projectdirs/cosmo/staging/'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-90prime-dr9.fits.gz'
plot_dir = '/global/u2/r/rongpu/temp/90prime_junk/brick_plots'

# Load brick list
bricks = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/survey-bricks.fits.gz')

# Load the list of CCDs left to process
ccd_columns = ['image_filename', 'image_hdu', 'ccdname', 'expnum', 'filter', 'mjd_obs', 'ra', 'dec', 'ccdraoff', 'ccddecoff', 'ccd_cuts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))

with open("/global/u2/r/rongpu/temp/list_of_ccds.txt", "r") as f:
    ccd_info_list = f.read().splitlines()

# brick_list = []
# for index in range(len(ccd_info_list)):

def find_overlapping_bricks(ccd_info):

    expnum, ccdname = ccd_info.split()
    expnum = int(expnum)
    print(expnum, ccdname)

    mask = (ccd['expnum']==expnum) & (ccd['ccdname']==ccdname)
    ccd_index = np.where(mask)[0][0]

    # Find bricks that cover the CCDs in this exposure
    search_radius = 1.414*((pixel_size*4079/2)+(0.26*3600/2)) # CCD size + brick size
    idx1, idx2, d2d, d_ra, d_dec = search_around(ccd['ra'][[ccd_index]], ccd['dec'][[ccd_index]], bricks['RA'], bricks['DEC'], search_radius=search_radius, verbose=False)
    mask = np.abs(d_ra)<(pixel_size*4079/2+0.26*3600/2)
    mask &= np.abs(d_dec)<(pixel_size*4079/2+0.26*3600/2)
    idx1, idx2 = idx1[mask], idx2[mask]

    brick_idx = idx2
    brick_list = list(bricks['BRICKNAME'][brick_idx])

    image_path = os.path.join(image_dir, ccd['image_filename'][ccd_index]).strip()
    hdu_index = ccd['image_hdu'][ccd_index]
    with fits.open(image_path) as hdulist:
        w = wcs.WCS(hdulist[hdu_index].header)
        naxis1 = hdulist[hdu_index].header['NAXIS1']
        naxis2 = hdulist[hdu_index].header['NAXIS2']

    # here x and y are not numpy indices, but they are the x-y values in the DECam CCD schematics
    pix_x_grid, pix_y_grid = np.meshgrid(np.arange(naxis1), np.arange(naxis2))
    pix_x, pix_y = pix_x_grid.flatten(), pix_y_grid.flatten()
    pix_ra, pix_dec = w.wcs_pix2world(pix_x, pix_y, 0)
    pix_ra = pix_ra + ccd['ccdraoff'][ccd_index] / 3600.
    pix_dec = pix_dec + ccd['ccddecoff'][ccd_index] / 3600.

    points = w.wcs_pix2world([[0,0], [naxis1, naxis2], [0, naxis2], [naxis1, 0]], 0)
    plt.figure(figsize=(6, 6))
    plt.plot(pix_ra[::500], pix_dec[::500], '.', ms=0.5)
    for ii in brick_idx:
        ra_plot = [bricks['RA1'][ii], bricks['RA1'][ii], bricks['RA2'][ii], bricks['RA2'][ii], bricks['RA1'][ii]]
        dec_plot = [bricks['DEC1'][ii], bricks['DEC2'][ii], bricks['DEC2'][ii], bricks['DEC1'][ii], bricks['DEC1'][ii]]
        plt.plot(ra_plot, dec_plot, alpha=0.5)
    plt.plot(points[:, 0], points[:, 1], '.', ms=5, color='r')
    plt.xlabel('RA')
    plt.ylabel('DEC')
    plt.title('{}-{}  {} overlapping bricks'.format(expnum, ccdname, len(brick_idx)))
    # plt.tight_layout()
    plt.savefig(os.path.join(plot_dir, '{}_{}.png'.format(expnum, ccdname)))
    plt.close()

    return brick_list


def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(find_overlapping_bricks, ccd_info_list)

    print('All done!!!!!!!!!!!!!!!')
    brick_list = []
    for tmp in res:
        brick_list += tmp

    brick_list_unique = np.unique(brick_list)
    
    with open("list_of_bricks.txt", "w") as f:
        for tmp in brick_list_unique:
            print(tmp)
            f.write(tmp+'\n')


if __name__=="__main__":
    main()
