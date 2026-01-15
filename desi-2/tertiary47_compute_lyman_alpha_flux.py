from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio


dirname_daily = '/global/cfs/cdirs/desi/spectro/redux/darknight/tiles/cumulative/83577/20250430'

version = 'laelbg_templates_2.89_3.19_archetype_galaxy_only'
dirname_rr = os.path.join('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary47', version)

fns = glob.glob(os.path.join(dirname_rr, 'redrock*.fits'))
# print(len(fns))
cat = []
for fn in fns:
    cat = Table(fitsio.read(fn, ext='REDSHIFTS'))

    coadd_fn = os.path.join(dirname_daily, os.path.basename(fn).replace('redrock-', 'coadd-'))
    wave = fitsio.read(coadd_fn, ext='B_WAVELENGTH')
    flux = fitsio.read(coadd_fn, ext='B_FLUX')

    cat['lya_flux'] = -99.  # in 10^-17 erg/cm^2/s
    for index in range(len(cat)):
        wave_min, wave_max = (1215.67-4.5) * (1+cat['Z'][index]), (1215.67+4.5) * (1+cat['Z'][index])
        mask = (wave>wave_min) & (wave<wave_max)
        cat['lya_flux'][index] = np.sum(flux[index][mask] * (wave[1]-wave[0]))

    cat = cat[['TARGETID', 'lya_flux']]
    cat.write(fn.replace('redrock-', 'lyaflux-'))
