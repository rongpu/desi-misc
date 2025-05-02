from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
from multiprocessing import Pool


stack_type = 'same_night'
# stack_type = 'different nights 3 months'


if stack_type == 'same_night':
    spec_fn = '/pscratch/sd/r/rongpu/tmp/sky_spectra/sky_spectra_same_night_20240212.fits'
elif stack_type == 'different nights 3 months':
    spec_fn = '/pscratch/sd/r/rongpu/tmp/sky_spectra/sky_spectra_different_nights_202401_202403.fits'

# spec = Table(fitsio.read(spec_fn))

cat = Table(fitsio.read(spec_fn, columns=['FIBER', 'TILEID']))

# Only keep unique fibers
_, idx = np.unique(cat['FIBER'], return_index=True)
spec = Table(fitsio.read(spec_fn, rows=idx))
assert len(spec)==500

fn_original = '/global/cfs/cdirs/desi/spectro/redux/loa/tiles/cumulative/1263/20240212/coadd-0-1263-thru20240212.fits'
fn_new = '/pscratch/sd/r/rongpu/tmp/sky_spectra/stacked_sky/test/coadd-0-1263-thru20240212-new-sky.fits'

fits = fitsio.FITS(fn_original)

# wave = {}
# for camera in ['B', 'R', 'Z']:
#     wave[camera] = fitsio.read(fn_original, ext=camera+'_WAVELENGTH')

with fitsio.FITS(fn_new, 'rw') as f:
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
        elif '_WAVELENGTH' not in extname and extname!='':
            data = np.array(spec[extname])
        f.write(data, extname=extname, header=header)

fits.close()


