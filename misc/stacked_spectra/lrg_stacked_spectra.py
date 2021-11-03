#!/usr/bin/env python

# import matplotlib
# matplotlib.use('Agg')
import sys
import os
from glob import glob
import numpy as np
import fitsio
import astropy.io.fits as fits
import fitsio
from desitarget.targetmask import desi_mask
from desitarget.geomask import match
# from raichoorlib import get_line_list
# import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.ticker import MultipleLocator

import time
from astropy.table import Table, vstack, hstack, join


time_start = time.time()

rfwmin, rfwmax, delta = 1800., 7017., 0.1
rfws = np.round(np.arange(rfwmin, rfwmax + delta, delta), 2)

d = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/tmp/lrgs_for_stacking.fits'))

# ###################################
# np.random.seed(718751)
# idx = np.random.choice(len(d), size=1000, replace=False)
# d = d[idx]
# ###################################

n = np.zeros(len(rfws), dtype=int)
fl, iv = np.zeros(len(rfws)), np.zeros(len(rfws))
#
tileids = np.unique(d['TILEID'])

for index, tileid in enumerate(tileids):

    # targetid and zspec
    sel0 = d['TILEID']==tileid
    thrunight = d['thrunight'][sel0][0]
    petals = np.unique(d['PETAL_LOC'][sel0])

    print(index, '/', len(tileids), sel0.sum())

    for petal in petals:
        sel = (d['TILEID']==tileid) & (d['PETAL_LOC']==petal)
        tids = d['TARGETID'][sel]
        zs = d['Z'][sel]

        fn = '/global/cfs/cdirs/desi/spectro/redux/everest/tiles/cumulative/{}/{}/coadd-{}-{}-thru{}.fits'.format(tileid, thrunight, petal, tileid, thrunight)

        h = fits.open(fn)
        ii, iisp = match(tids, h['FIBERMAP'].data['TARGETID'])
        # looping on the cameras
        for camera in ['B', 'R', 'Z']:
            # reading
            wsp = h['{}_WAVELENGTH'.format(camera)].data
            flsp = h['{}_FLUX'.format(camera)].data
            ivsp = h['{}_IVAR'.format(camera)].data
            # looping through each spectrum...
            # interpolating to rest-frame wavelengths
            for i, isp in zip(ii, iisp):
                rfflsp = np.interp(rfws, wsp / (1 + zs[i]), flsp[isp], left=0, right=0)
                rfivsp = np.interp(rfws, wsp / (1 + zs[i]), ivsp[isp], left=0, right=0)
                ws = rfivsp
                fl += rfflsp * ws
                iv += ws
                n[ws > 0] += 1
#
sel = iv > 0
fl[sel] /= iv[sel]
#
cols = []
cols += [fits.Column(name='WAVELENGTH', format='E', array=rfws)]
cols += [fits.Column(name='FLUX', format='E', array=fl)]
cols += [fits.Column(name='NSPEC', format='K', array=n)]
h = fits.BinTableHDU.from_columns(fits.ColDefs(cols))
h.writeto('/global/cfs/cdirs/desi/users/rongpu/tmp/lrg_stacked_spectra.fits', overwrite=True)

print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
