# salloc -N 1 -C haswell -q interactive -t 04:00:00
# shifter --image docker:legacysurvey/legacypipe:DR9.5.5 bash

# Next: 1. require a minimum of 20 sources to qualify
#       2. change threshold2 to 3.0

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from astropy import units as u
from astropy.coordinates import SkyCoord

from multiprocessing import Pool
from pathlib import Path

##################################################################################################################

n_node = 4
task_id = 3

cleanup = True

plot_q = True
plot_dir = '/global/cscratch1/sd/rongpu/temp/double_vision_plots_round_2'

image_dir = '/global/project/projectdirs/cosmo/staging'
tmp_dir = '/global/cscratch1/sd/rongpu/temp/double_vision'

# Include spaces for DR9 CCDs, e.g. 'N9 '
ccds = ['S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'S8', 'S9', 'S10', 'S11', 'S12', 'S13', 'S14', 'S15', 
        'S16', 'S17', 'S18', 'S19', 'S20', 'S21', 'S22', 'S23', 'S24', 'S25', 'S26', 'S27', 'S28', 'S29', 
        'S31', 'N1', 'N2', 'N3', 'N4', 'N5', 'N6', 'N7', 'N8', 'N9', 'N11', 'N12', 'N14', 'N16', 'N17', 
        'N18', 'N19', 'N20', 'N21', 'N22', 'N23', 'N24', 'N25', 'N26', 'N27', 'N28', 'N29', 'N31']
# ccds_more = []
n_ccd = 3
n_ccd_threshold = 2 # require at least this many CCDs to meet the cut for visual inspection
search_radius = 30. # arcsec

threshold1 = 2.
threshold2 = 4.

n_processes = 32

image_vrange = {'g':15, 'r':20, 'z':50}

##################################################################################################################

surveyccd_path_dr9 = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'

ccd_columns = ['expnum', 'image_filename', 'ccdname', 'filter', 'ccd_cuts']
ccd = Table(fitsio.read(surveyccd_path_dr9, columns=ccd_columns))
print(len(ccd))

def search_around(ra1, dec1, ra2, dec2, search_radius=1., verbose=True):
    '''
    Using the astropy.coordinates.search_around_sky module to find all pairs within
    some search radius.

    Inputs: 
    RA and Dec of two catalogs;
    search_radius (arcsec);


    Outputs: 
        idx1, idx2: indices of matched objects in the two catalogs;
        d2d: angular distances (arcsec);
        d_ra, d_dec: the differences in RA and Dec (arcsec); 
    '''
    
    # protect the global variables from being changed by np.sort
    ra1, dec1, ra2, dec2 = map(np.copy, [ra1, dec1, ra2, dec2])
    
    # Matching catalogs
    sky1 = SkyCoord(ra1*u.degree,dec1*u.degree, frame='icrs')
    sky2 = SkyCoord(ra2*u.degree,dec2*u.degree, frame='icrs')
    idx1, idx2, d2d, d3d = sky2.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
    if verbose:
        print('%d nearby objects'%len(idx1))
    
    # convert distances to numpy array in arcsec
    d2d = np.array(d2d.to(u.arcsec))


    d_ra = (ra2[idx2]-ra1[idx1])*3600.    # in arcsec
    d_dec = (dec2[idx2]-dec1[idx1])*3600. # in arcsec

    ##### Convert d_ra to actual arcsecs #####
    mask = d_ra > 180*3600
    d_ra[mask] = d_ra[mask] - 360.*3600
    mask = d_ra < -180*3600
    d_ra[mask] = d_ra[mask] + 360.*3600
    d_ra = d_ra * np.cos(dec1[idx1]/180*np.pi)
    ##########################################

    return idx1, idx2, d2d, d_ra, d_dec


def find_candidates(expnum):

    ccd_threshold_counter = 0
        
    idx = np.where(ccd['expnum']==expnum)[0]
    idx = idx[np.in1d(ccd['ccdname'][idx], ccds)]
        
    if len(idx)<n_ccd_threshold:
        print('Not enough CCDs; skip')
        return None
    elif len(idx)>n_ccd:
        np.random.seed(expnum)
        ccds1 = np.random.choice(np.unique(ccd['ccdname'][idx]), size=n_ccd, replace=False)
        idx = idx[np.in1d(ccd['ccdname'][idx], ccds1)]

    for ccd_index in idx:

        meet_criteria = False

        ccdname = ccd['ccdname'][ccd_index].strip()
        band = ccd['filter'][ccd_index]

        image_path = os.path.join(image_dir, ccd['image_filename'][ccd_index].strip())
        downsized_image_path = os.path.join(tmp_dir, 'image_decam_{}_{}_{}_small.fits'.format(expnum, ccdname, band))

        print(expnum, ccdname, band)

        ############################ run image2xy ############################

        cat_path = downsized_image_path.replace('small.fits', 'small.xy.fits')

        img = fitsio.read(image_path, ext=ccdname)

        # Apply masks for DR9 images
        ood_path = image_path.replace('_ooi_', '_ood_')
        ood = fitsio.read(ood_path, ext=ccdname)
        img[ood!=0] = np.nan
            
        # downsize image
        binsize = 8
        trim_size_x = img.shape[1] % binsize
        trim_size_y = img.shape[0] % binsize
        img = img[:(img.shape[0]-trim_size_y), :(img.shape[1]-trim_size_x)]
        # to ignore NAN values, use np.nanmean
        # img = np.mean(np.mean(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)
        img = np.nanmean(np.nanmean(img.reshape((img.shape[0]//binsize, binsize, img.shape[1]//binsize,-1)), axis=3), axis=1)
        
        mask = ~np.isfinite(img)
        img[mask] = 0

        fitsio.write(downsized_image_path, img)
        
        if os.path.isfile(cat_path):
            print(cat_path+' already exists! skip image2xy')
        else:
            os.system('image2xy '+downsized_image_path)

        if cleanup:
            os.remove(downsized_image_path)

        ####################################################################

        try:
            cat = Table(fitsio.read(cat_path))
        except:
            print('Error loading catalog!!!')
            return None
        cat['RA'], cat['DEC'] = cat['Y']*0.262/3600, cat['X']*0.262/3600

        idx1, idx2, d2d, d_ra, d_dec = search_around(cat['RA'], cat['DEC'], cat['RA'], cat['DEC'], search_radius=search_radius, verbose=False)
        mask = (d_ra**2+d_dec**2)>0.5
        d_ra = d_ra[mask]
        d_dec = d_dec[mask]

        xbins, ybins = np.linspace(-search_radius, search_radius, 50), np.linspace(-search_radius, search_radius, 50)
        counts, _, _, = np.histogram2d(d_ra, d_dec, bins=[xbins, ybins])
        mask = counts>threshold1 * (np.percentile(counts.flatten(), 99)) # candidate criteria
        if np.sum(mask)>=2:
            meet_criteria = True
            if plot_q:
                fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
                counts, _, _, im, = ax[1].hist2d(d_ra, d_dec, bins=[xbins, ybins])
                ax[0].hist(counts.flatten(), 100)
                ax[0].axvline(np.percentile(counts.flatten(), 99), color='r', lw=1)
                ax[0].axvline(threshold1*np.percentile(counts.flatten(), 99), color='r', lw=1)
                # ax[1].axhline(0, color='white', lw=1)
                # ax[1].axvline(0, color='white', lw=1)
                ax[1].set_xlabel('d_ra')
                ax[1].set_ylabel('d_dec')
                plt.colorbar(im)
                plt.savefig(os.path.join(plot_dir, 'image_decam_{}_{}_{}_hist2d_1.jpeg'.format(expnum, ccdname, band)))
                plt.close()

        # a different bin size
        xbins, ybins = np.linspace(-search_radius, search_radius, 120), np.linspace(-search_radius, search_radius, 120)
        counts, _, _, = np.histogram2d(d_ra, d_dec, bins=[xbins, ybins])
        mask = counts > threshold2 * (np.percentile(counts.flatten(), 99)) # candidate criteria
        if np.sum(mask)>=2:
            meet_criteria = True
            if plot_q:
                fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
                counts, _, _, im, = ax[1].hist2d(d_ra, d_dec, bins=[xbins, ybins])
                ax[0].hist(counts.flatten(), 100)
                ax[0].axvline(np.percentile(counts.flatten(), 99), color='r', lw=1)
                ax[0].axvline(threshold2*np.percentile(counts.flatten(), 99), color='r', lw=1)
                # ax[1].axhline(0, color='white', lw=1)
                # ax[1].axvline(0, color='white', lw=1)
                ax[1].set_xlabel('d_ra')
                ax[1].set_ylabel('d_dec')
                plt.colorbar(im)
                plt.savefig(os.path.join(plot_dir, 'image_decam_{}_{}_{}_hist2d_2.jpeg'.format(expnum, ccdname, band)))
                plt.close()

        if meet_criteria:
            ccd_threshold_counter += 1

        if meet_criteria and plot_q:
            vrange = image_vrange[band]
            # naive sky estimation
            mask = (img<np.percentile(img.flatten(), 95))
            median_sky = np.median(img[mask].flatten())

            plt.figure(figsize=(12, 6))
            plt.imshow((img-median_sky).T, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, 'image_decam_{}_{}_{}_small.jpeg'.format(expnum, ccdname, band)))
            plt.close()

            plt.figure(figsize=(12, 6))
            plt.imshow((img-median_sky).T, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
            plt.plot(cat['Y']-1, cat['X']-1, 'bx', ms=12, alpha=0.5)
            plt.tight_layout()
            plt.savefig(os.path.join(plot_dir, 'image_decam_{}_{}_{}_sources.jpeg'.format(expnum, ccdname, band)))
            plt.close()

        # if cleanup and (meet_criteria==False):
        if cleanup:
            os.remove(cat_path)

        if ccd_threshold_counter>=n_ccd_threshold:
            break

    if ccd_threshold_counter>=n_ccd_threshold:
        print('Found a candidate: expnum {}'.format(expnum))
        if expnum in known_bad_expnum_list:
            Path(os.path.join(tmp_dir, 'known_decam_{}'.format(expnum))).touch()
        else:
            Path(os.path.join(tmp_dir, 'decam_{}'.format(expnum))).touch()

    return None

##################################

# expnum_list = [425619, 548188, 548218, 548362, 548400, 747320, 548195]
# for expnum in expnum_list:
#     find_candidates(expnum)

# np.random.seed(623)
# expnum_list = np.concatenate([expnum_list, np.random.choice(ccd['expnum'], 100, replace=False)])

# np.random.seed(623)
# expnum_list = np.random.choice(ccd['expnum'], 100, replace=False)
# # expnum_list = np.unique(expnum_list)

# for expnum in expnum_list:
#     find_candidates(expnum)

expnum_list = np.unique(ccd['expnum'])

# mask = ccd['ccd_cuts']==0
# expnum_list = np.unique(ccd['expnum'][mask])
# print(len(expnum_list))

# ###############
# np.random.seed(623)
# expnum_list = np.random.choice(expnum_list, 640, replace=False)
# ###############

expnum_list_split = np.array_split(expnum_list, n_node)
expnum_list = np.sort(expnum_list_split[task_id])

with open('decam-bad_expid_20200814.txt') as f:
    texts = list(map(str.strip, f.readlines()))
known_bad_expnum_list = []
for text in texts:
    if len(text)>0 and (text[0]!='#'):
        known_bad_expnum_list.append(int(text.replace('-', ' ').split()[0]))

start = time.time()
with Pool(processes=n_processes) as pool:
    res = pool.map(find_candidates, expnum_list)
end = time.time()
print('Took {:.1f} seconds'.format(end - start))

# def main():

#     start = time.time()

#     with Pool(processes=n_processes) as pool:
#         res = pool.map(find_candidates, expnum_list)

#     # make_plots(expnum_list[0])

#     end = time.time()
#     print('Took {:.1f} seconds'.format(end - start))

#     print('Done!!!!!!!!!!!!!!!!!!!!!')

# if __name__=="__main__":
#     main()

