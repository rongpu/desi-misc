from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio


fns = glob.glob('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/daily/rr_lae_2.25_3.4_nmf/redrock*.fits')
print(len(fns))

cat = []
for fn in fns:

    coadd_fn = fn.replace('redrock-', 'coadd-').replace('/rr_lae_2.25_3.4_nmf/', '/')

    cat = Table(fitsio.read(fn, ext='REDSHIFTS'))

    wave = fitsio.read(coadd_fn, ext='B_WAVELENGTH')
    flux = fitsio.read(coadd_fn, ext='B_FLUX')
    flux_var = 1/(fitsio.read(coadd_fn, ext='B_IVAR'))

    cat['lya_flux'] = -99.  # unit: 10^-17 erg/cm^2/s
    cat['lya_flux_err'] = -99.  # unit: 10^-17 erg/cm^2/s
    for index in range(len(cat)):
        wave_min, wave_max = (1215.67-4.5) * (1+cat['Z'][index]), (1215.67+4.5) * (1+cat['Z'][index])
        mask = (wave>wave_min) & (wave<wave_max)
        cat['lya_flux'][index] = np.sum(flux[index][mask]) * (wave[1]-wave[0])
        cat['lya_flux_err'][index] = np.sqrt(np.sum(flux_var[index][mask])) * (wave[1]-wave[0])

    cat = cat[['TARGETID', 'lya_flux', 'lya_flux_err']]
    lya_path = fn.replace('redrock-', 'lyaflux-')
    cat.write(lya_path, overwrite=True)
