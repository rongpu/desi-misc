from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

bricks = Table(fitsio.read('/global/cfs/cdirs/cosmo/work/users/rongpu/survey-bricks-10x10-decam-deep-fields-cosmos.fits'))
print(len(bricks))

np.random.seed(98175)
idx = np.random.choice(len(bricks), size=30, replace=False)
bricks = bricks[idx]

bricknames = np.unique(bricks['BRICKNAME'])
brickname_parent_list = [brickname[:-4] for brickname in bricknames]
brickname_parent_list = np.unique(brickname_parent_list)
print(len(brickname_parent_list))


def create_image(data, cmap='gray', dpi=80, vmin=None, vmax=None, origin=None, norm=None):
    '''
    Create an image with exactly the same pixel dimension as the data.
    Example:
        x, y = np.arange(0, 10), np.arange(0, 10)
        xx, yy = np.meshgrid(x, y)
        img = np.array((xx + yy)%2==0, dtype=int)
        ax = create_image(img)
        plt.savefig('img.png')
        plt.close()
    '''
    xpixels, ypixels = data.shape[0], data.shape[1]
    figsize = ypixels / dpi, xpixels / dpi
    fig = plt.figure(figsize=figsize, dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(data, cmap=cmap, interpolation='none', vmin=vmin, vmax=vmax, origin=origin, norm=norm)
    plt.axis('off')
    return ax


vranges = {'g':0.015, 'r':0.03, 'i':0.04, 'z':0.05, 'W1':15, 'W2':40}

for brickname_parent in brickname_parent_list:

    for band in ['g', 'r', 'i', 'z']:
        vrange = vranges[band]
        vrange = vrange/2

        image_mosaic = np.zeros((3440, 3440))
        for i in range(0, 10):
            for j in range(0, 10):
                brickname = brickname_parent+'.{}.{}'.format(i, j)
                coadd_dir = '/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/coadd/{}/{}'.format(brickname[:3], brickname)
                metrics_dir = '/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/metrics/{}'.format(brickname[:3])
                image_path = os.path.join(coadd_dir, 'legacysurvey-{}-image-{}.fits.fz'.format(brickname, band))
                if os.path.isfile(image_path):
                    image = fitsio.read(image_path)
                else:
                    # print('Missing:', brickname)
                    continue
                image_mosaic[344*j:344*(j+1), 344*(9-i):344*(9-i+1)] = image[76:420, 76:420]
        # image_mosaic = image_downsize(image_mosaic, binsize=2)
        create_image(image_mosaic, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
        plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/coadd_plots/dr10_deep/deep_{}_{}_image.png'.format(brickname_parent, band))
        plt.close()

        model_mosaic = np.zeros((3440, 3440))
        for i in range(0, 10):
            for j in range(0, 10):
                brickname = brickname_parent+'.{}.{}'.format(i, j)
                coadd_dir = '/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/coadd/{}/{}'.format(brickname[:3], brickname)
                metrics_dir = '/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/metrics/{}'.format(brickname[:3])
                model_path = os.path.join(coadd_dir, 'legacysurvey-{}-model-{}.fits.fz'.format(brickname, band))
                if os.path.isfile(model_path):
                    model = fitsio.read(model_path)
                else:
                    # print('Missing:', brickname)
                    continue
                model_mosaic[344*j:344*(j+1), 344*(9-i):344*(9-i+1)] = model[76:420, 76:420]
        # model_mosaic = model_downsize(model_mosaic, binsize=2)
        create_image(model_mosaic, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
        plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/coadd_plots/dr10_deep/deep_{}_{}_model.png'.format(brickname_parent, band))
        plt.close()

        resid_mosiac = image_mosaic - model_mosaic
        create_image(resid_mosiac, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
        plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/coadd_plots/dr10_deep/deep_{}_{}_resid.png'.format(brickname_parent, band))
        plt.close()


    # ######################################## DR9 ############################################

    # coadd_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/coadd/{}/{}'.format(brickname_parent[:3], brickname_parent)
    # metrics_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/metrics/{}'.format(brickname_parent[:3])

    # for band in ['g', 'r', 'z']:
    #     vrange = vranges[band]
    #     vrange = vrange/2

    #     image_path = os.path.join(coadd_dir, 'legacysurvey-{}-image-{}.fits.fz'.format(brickname_parent, band))
    #     image = fitsio.read(image_path)
    #     image = image[80:3520, 80:3520]
    #     # image = image_downsize(image, binsize=2)
    #     create_image(image, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
    #     plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/coadd_plots/dr9/dr9_{}_{}_image.png'.format(brickname_parent, band))
    #     plt.close()

    #     model_path = os.path.join(coadd_dir, 'legacysurvey-{}-model-{}.fits.fz'.format(brickname_parent, band))
    #     model = fitsio.read(model_path)
    #     model = model[80:3520, 80:3520]
    #     # model = model_downsize(model, binsize=2)
    #     create_image(model, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
    #     plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/coadd_plots/dr9/dr9_{}_{}_model.png'.format(brickname_parent, band))
    #     plt.close()

    #     resid = image - model
    #     create_image(resid, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
    #     plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/coadd_plots/dr9/dr9_{}_{}_resid.png'.format(brickname_parent, band))
    #     plt.close()


    # ######################################## DR9 ############################################

    # coadd_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9.1.1/coadd/{}/{}'.format(brickname_parent[:3], brickname_parent)
    # metrics_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9.1.1/metrics/{}'.format(brickname_parent[:3])

    # for band in ['g', 'r', 'z']:
    #     vrange = vranges[band]
    #     vrange = vrange/2

    #     image_path = os.path.join(coadd_dir, 'legacysurvey-{}-image-{}.fits.fz'.format(brickname_parent, band))
    #     image = fitsio.read(image_path)
    #     image = image[80:3520, 80:3520]
    #     # image = image_downsize(image, binsize=2)
    #     create_image(image, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
    #     plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/coadd_plots/dr9.1.1/dr9.1.1_{}_{}_image.png'.format(brickname_parent, band))
    #     plt.close()

    #     model_path = os.path.join(coadd_dir, 'legacysurvey-{}-model-{}.fits.fz'.format(brickname_parent, band))
    #     model = fitsio.read(model_path)
    #     model = model[80:3520, 80:3520]
    #     # model = model_downsize(model, binsize=2)
    #     create_image(model, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
    #     plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/coadd_plots/dr9.1.1/dr9.1.1_{}_{}_model.png'.format(brickname_parent, band))
    #     plt.close()

    #     resid = image - model
    #     create_image(resid, cmap='seismic', vmin=-vrange, vmax=vrange, origin='lower')
    #     plt.savefig('/global/cfs/cdirs/cosmo/www/temp/rongpu/dr10dev/deep_fields/coadd_plots/dr9.1.1/dr9.1.1_{}_{}_resid.png'.format(brickname_parent, band))
    #     plt.close()