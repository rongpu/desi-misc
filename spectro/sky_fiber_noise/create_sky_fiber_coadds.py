from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
from multiprocessing import Pool


b_only = False
br_only = True

# stack_type = 'same night'
# # stack_type = 'different nights 3 months'

# n_stack = 1


# for stack_type in ['same night', 'different nights 3 months', 'same night shuffle fibers']:
for stack_type in ['same night shuffle fibers']:

    fn_original = '/global/cfs/cdirs/desi/spectro/redux/loa/tiles/cumulative/1263/20240212/coadd-0-1263-thru20240212.fits'

    if (stack_type == 'same night') or (stack_type == 'same night shuffle fibers'):
        spec_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/spectro/sky_spectra/sky_spectra_same_night_20240212.fits'
    elif stack_type == 'different nights 3 months':
        spec_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/spectro/sky_spectra/sky_spectra_different_nights_202401_202403.fits'

    # spec = Table(fitsio.read(spec_fn))

    specall = Table(fitsio.read(spec_fn))
    print(len(specall))

    if stack_type == 'same night shuffle fibers':
        # shuffle
        np.random.seed(111)
        idx = np.random.choice(len(specall), size=len(specall), replace=False)
        specall = specall[idx]
        print(len(specall))

        specall['FIBER_ORIG'] = specall['FIBER'].copy()
        fibers = np.repeat(np.arange(500), 32)
        specall = specall[:len(fibers)]
        specall['FIBER'] = fibers
        print(len(specall), len(np.unique(specall['FIBER'])))

        specall['keep'] = False
        for fiber in np.unique(specall['FIBER']):
            idx = np.where(specall['FIBER']==fiber)[0]
            if len(np.unique(specall['FIBER_ORIG'][idx]))<len(idx):
                _, idx1 = np.unique(specall['FIBER_ORIG'][idx], return_index=True)
                idx = idx[idx1]
            specall['keep'][idx] = True
        specall = specall[specall['keep']]
        print(len(specall), len(np.unique(specall['FIBER'])))

        fibercount = Table()
        fibercount['FIBER'], fibercount['count'] = np.unique(specall['FIBER'], return_counts=True)

    fibers = np.unique(specall['FIBER'])
    print(len(fibers))

    fibercount = Table()
    fibercount['FIBER'], fibercount['count'] = np.unique(specall['FIBER'], return_counts=True)
    print('minimum number of spectra per fiber: ', fibercount['count'].min())
    print('maximum number of spectra per fiber: ', fibercount['count'].max())

    for n_stack in [1, 4, 8, 32]:
        if stack_type == 'same night shuffle fibers' and n_stack==1:
            continue

        if stack_type == 'same night':
            fn_output = '/global/cfs/cdirs/desicollab/users/rongpu/data/spectro/sky_spectra/stacked_sky/coadd-same_night_20240212-stack_{}.fits'.format(n_stack)
        elif stack_type == 'different nights 3 months':
            fn_output = '/global/cfs/cdirs/desicollab/users/rongpu/data/spectro/sky_spectra/stacked_sky/coadd-different_nights_202401_202403-stack_{}.fits'.format(n_stack)
        elif stack_type == 'same night shuffle fibers':
            fn_output = '/global/cfs/cdirs/desicollab/users/rongpu/data/spectro/sky_spectra/stacked_sky/coadd-same_night_shuffle_fibers_20240212-stack_{}.fits'.format(n_stack)

        if b_only:
            fn_output = fn_output.replace('.fits', '-b_only.fits')
        elif br_only:
            fn_output = fn_output.replace('.fits', '-br_only.fits')

        if os.path.isfile(fn_output):
            print('file already exists; skip')
            print(fn_output)
            continue

        spec_stack = []
        for fiber in fibers:
            spec = Table()
            spec['FIBER'] = [fiber]
            idx = np.where(specall['FIBER']==fiber)[0]
            spec1 = specall[idx]
            for camera in ['B', 'R', 'Z']:
                ffivar = np.sum(spec1[camera+'_IVAR'][:n_stack] * spec1[camera+'_MSK'][:n_stack], axis=0)
                ff = np.sum(spec1[camera+'_FLUX'][:n_stack] * spec1[camera+'_IVAR'][:n_stack] * spec1[camera+'_MSK'][:n_stack], axis=0) / ffivar
                ff[np.isnan(ff)] = 0.
                msk = np.sum(np.array(spec1[camera+'_MSK'][:n_stack]), axis=0)>0
                spec[camera+'_FLUX'] = [ff]
                spec[camera+'_IVAR'] = [ffivar]
                spec[camera+'_MSK'] = [msk]
                spec[camera+'_RESOLUTION'] = [spec1[camera+'_RESOLUTION'][0]]
            spec_stack.append(spec)
        spec = vstack(spec_stack, join_type='exact')

        fits = fitsio.FITS(fn_original)

        with fitsio.FITS(fn_output, 'rw') as f:
            for index in range(len(fits)):
                extname = fits[index].get_extname()
                data, header = fitsio.read(fn_original, ext=index, header=True)
                if extname in ['FIBERMAP', 'EXP_FIBERMAP', 'SCORES']:
                    for col in data.dtype.names:
                        data[col] = data[col][0]
                    data['TARGETID'] = np.arange(500)
                    if extname=='FIBERMAP':
                        data['FIBER'] = np.array(spec['FIBER'])
                elif '_MASK' in extname:
                    data = np.array(~spec[extname[0]+'_MSK']).astype(np.int32) * 2**1  # ALLBADPIX
                elif extname[1:] in ['_FLUX', '_IVAR', '_MSK', '_RESOLUTION']:
                    data = np.array(spec[extname])
                if b_only and (extname.startswith('R_') or extname.startswith('Z_')):  # skip R and Z cameras
                    continue
                if br_only and extname.startswith('Z_'):  # skip Z camera
                    continue
                f.write(data, extname=extname, header=header)

        fits.close()


