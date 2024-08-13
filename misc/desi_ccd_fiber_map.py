from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from desispec.io import read_xytraceset


fiber_pix = 2  # width of the fiber trace is fiber_pix*2+1 pixels (by default fiber_pix=2 -> 5 pixels wide)

psf_dir = '/global/cfs/cdirs/desi/spectro/redux/jura/exposures/20230423/00177541'
output_dir = '/global/cfs/cdirs/desicollab/users/rongpu/useful/desi_ccd_fiber_map_20230423/'

# expid = psf_dir.split('/')[-1]

fm_dict = {}

for petal in range(10):
    print(petal)
    for camera in ['b', 'r', 'z']:

        cp = camera+str(petal)
        psf_fn = os.path.join(psf_dir, 'psf-{}-00177541.fits'.format(cp))
        tset = read_xytraceset(psf_fn)

        fibers = np.arange(500)

        nw = 50
        wave = np.linspace(tset.wavemin, tset.wavemax, nw)

        if camera=='b':
            fibermap = np.ones((4096, 4096), dtype=np.int32) * -1
        else:
            fibermap = np.ones((4128, 4114), dtype=np.int32) * -1
        # CCD first index: wavelength direction; second index: fiber index direction

        for fiber in fibers:

            actual_fiber = fiber + petal * 500

            x = tset.x_vs_wave(fiber, wave)
            y = tset.y_vs_wave(fiber, wave)

            xmin, xmax = 0, fibermap.shape[1]-1
            ymin, ymax = 0, fibermap.shape[0]-1

            x_raw = tset.x_vs_wave(fiber, wave)
            y_raw = tset.y_vs_wave(fiber, wave)

            y = np.arange(ymin, ymax+1)
            x = np.interp(y, y_raw, x_raw)

            x = np.round(x).astype(int)
            y = np.round(y).astype(int)
            mask = (x>=xmin) & (x<=xmax)
            mask &= (y>=ymin) & (y<=ymax)
            x = x[mask]
            y = y[mask]

            for ii in range(-fiber_pix, fiber_pix+1):
                assert np.sum(fibermap[y, x+ii]!=-1)==0, 'Fiber footprints overlap; try a smaller fiber_pix value'
                fibermap[y, x+ii] = actual_fiber

        fm_dict[cp] = fibermap.copy()

        output_path = os.path.join(output_dir, 'ccd_fibermap_{}.fits.gz'.format(cp))
        fits = fitsio.FITS(output_path, 'rw')
        fits.write(fibermap, extname=cp)
        fits.close()
