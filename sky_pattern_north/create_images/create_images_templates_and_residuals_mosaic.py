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

from scipy.ndimage.filters import gaussian_filter
from pathlib import Path
from scipy import stats
from multiprocessing import Pool

################################################################################

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'} 
plt.rcParams.update(params)

n_processes = 32

max_exposure = 50
vrange = 0.002
binsize = 2
pix_size = 0.262/3600*binsize
overwrite = True

################################################################################

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

ccdnamenumdict = {'CCD1': 1, 'CCD2': 2, 'CCD3': 3, 'CCD4': 4}
ccdnamenumdict_inv = {aa: bb for bb, aa in ccdnamenumdict.items()}

ccdnum_list = [1, 2, 3, 4]
ccd_ra = [-0.1554, -0.1554, 0.1554, 0.1554]
ccd_dec = [-0.1554, 0.1554, -0.1554, 0.1554]

################################################################################

image_dir = '/global/project/projectdirs/cosmo/staging'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/mosaic_ccd_blob_mask'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-mosaic-dr9.fits.gz'
template_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates_mosaic'

plot_dir = '/global/cfs/cdirs/desi/www/users/rongpu/plots/dr9dev/sky_pattern_north/sky_templates_v1/templates'

# skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8-mosaic.fits')
skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8-mosaic-subset.fits')
print('skyrun', len(skyrun))

mask = skyrun['ok']==True
skyrun = skyrun[mask]
print(len(skyrun))

sky_path_list = glob.glob(os.path.join(template_dir, '*.fits.fz'))
print(len(sky_path_list))

# ccd_columns = ['image_hdu', 'expnum', 'ccdname', 'ccdskycounts']
# ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))

def make_plots(sky_path):
    
    # # The file should be at least 0.5 hours old to ensure it's not being written
    # time_modified = os.path.getmtime(sky_path)
    # if (time.time() - time_modified)/3600 < 0.5:
    #     # continue
    #     return None

    run = int(sky_path[sky_path.rfind('_')+1:sky_path.rfind('.fits.fz')])
    
    # Get run info
    mask = skyrun['run']==run
    n_exposure = np.minimum(max_exposure, np.sum(mask))
    band = skyrun['filter'][mask][0]
    mjd_span = skyrun['mjd_obs'][mask].max() - skyrun['mjd_obs'][mask].min()

    if 'residual' in os.path.basename(sky_path):
        plot_type = 'residual'
    else:
        plot_type = 'template'

    plot_path = os.path.join(plot_dir, band, '{}_{}_sky_{}.png'.format(band, run, plot_type))

    if (overwrite==False) and os.path.isfile(plot_path):
        # print(plot_path, 'already exists!!!')
        return None

    if not os.path.exists(os.path.dirname(plot_path)):
        try: # in case another process is also creating the directory
            os.makedirs(os.path.dirname(plot_path))
        except:
            pass

    print(plot_path)
    Path(plot_path).touch()

    plt.figure(figsize=(8, 8))
    # for 90prime:
    # plt.figure(figsize=(11, 11))
    for ii, ccdnum in enumerate(ccdnum_list):

        # print(ii)
        ccdname = ccdnamenumdict_inv[ccdnum]

        try:
            sky = fits.getdata(sky_path, extname=ccdname)
        except:
            print(ccdname+' does not exist in image!')
            continue

        ################ downsize image ################

        # trim edges to enable downsizing
        # trimmed image size need to be multiples of binsize
        trim_size_x = sky.shape[1] % binsize
        trim_size_y = sky.shape[0] % binsize
        sky = sky[:(sky.shape[0]-trim_size_y), :(sky.shape[1]-trim_size_x)]

        if np.sum(~np.isfinite(sky))>0:
            print('{} bad pixels'.format(~np.isfinite(sky)))

        # to ignore NAN values, use np.nanmean
        sky = np.mean(np.mean(sky.reshape((sky.shape[0]//binsize, binsize, sky.shape[1]//binsize,-1)), axis=3), axis=1)

        ################################################

        if plot_type=='residual':
            # Gaussian filtering
            sky = gaussian_filter(sky, 2.5, mode='reflect')

        ysize, xsize = sky.shape
        ra, dec = ccd_ra[ii], ccd_dec[ii]

        fig = plt.imshow(sky.T, cmap='seismic', vmin=-vrange, vmax=vrange, 
                   extent=(ra-ysize*pix_size/2, ra+ysize*pix_size/2, dec+xsize*pix_size/2, dec-xsize*pix_size/2))

    text = 'run {}, {} band, {} exposures stacked, observed over {:.1f} days'.format(run, band, n_exposure, mjd_span)
    plt.title(text)
    plt.axis([0.3, -0.3, -0.3, 0.3])
    # for 90prime:
    # plt.axis([0.55, -0.55, -0.55, 0.55])
    plt.axis('off')
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    # plt.colorbar(fraction=0.04, pad=0.04)
    plt.tight_layout()
    plt.savefig(plot_path)
    plt.close()

def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(make_plots, sky_path_list)

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

