from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio


tile_list = [83579, 83580, 83581, 83582]
lastnight_list = [20260116, 20260115, 20260116, 20260116]

rr_dirname = '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/rr_lae_2.25_3.4_nmf'


fns = []
for tileid, lastnight in zip(tile_list, lastnight_list):
    print(tileid, lastnight)
    fns += glob.glob('/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/{}/{}/redrock-*thru{}.fits'.format(tileid, lastnight, lastnight))
fns.sort()
print(len(fns))


cat = []
for fn in fns:

    rr_fn = os.path.join(rr_dirname, os.path.basename(fn))
    coadd_fn = fn.replace('redrock-', 'coadd-')

    cat = Table(fitsio.read(rr_fn, ext='REDSHIFTS'))

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
    lya_path = os.path.join(rr_dirname, 'lyaflux', os.path.basename(rr_fn.replace('redrock-', 'lyaflux-')))
    cat.write(lya_path)
