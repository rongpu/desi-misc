from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
from multiprocessing import Pool


exposures = Table.read('/global/cfs/cdirs/desi/spectro/redux/daily/exposures-daily.csv')

fns = glob.glob('/global/cfs/cdirs/desi/spectro/redux/m4/exposures/*/*')
exp_list = []
for fn in fns:
    exp_list.append(int(fn.split('/')[-1]))
exp_list = np.sort(exp_list)
print(len(exp_list))

mask = np.in1d(exposures['EXPID'], exp_list)
exposures = exposures[mask]
print(len(exposures))

np.unique(exposures['FAFLAVOR'], return_counts=True)

mask = exposures['FAFLAVOR']=='maindark'
exposures = exposures[mask]
print(len(exposures))


camera = 'z'
np.random.seed(999)
idx = np.sort(np.random.choice(len(exposures), size=20, replace=False))

def process_exposure_petal(args):
    index, petal = args
    night = exposures['NIGHT'][index]
    expid = exposures['EXPID'][index]
    expid_str = str(expid).zfill(8)

    coadd_fn_new = f'/global/cfs/cdirs/desi/spectro/redux/m4/exposures/{night}/{expid_str}/cframe-{camera}{petal}-{expid_str}.fits.gz'
    if not os.path.isfile(coadd_fn_new):
        print('file does not exist:', coadd_fn_new)
        return

    coadd_fn_old = f'/global/cfs/cdirs/desi/spectro/redux/loa/exposures/{night}/{expid_str}/cframe-{camera}{petal}-{expid_str}.fits.gz'
    if not os.path.isfile(coadd_fn_old):
        coadd_fn_old = f'/global/cfs/cdirs/desi/spectro/redux/daily/exposures/{night}/{expid_str}/cframe-{camera}{petal}-{expid_str}.fits.gz'

    cat = Table(fitsio.read(coadd_fn_new, ext='FIBERMAP'))
    idx1 = np.where((cat['OBJTYPE'] == 'SKY') & (cat['FIBERSTATUS'] == 0))[0]
    n_plot = 2
    if len(idx1) > n_plot:
        idx1 = np.random.choice(idx1, size=n_plot, replace=False)

    wave = fitsio.read(coadd_fn_new, ext='WAVELENGTH')
    flux_new = fitsio.read(coadd_fn_new, ext='FLUX')
    msk_new = fitsio.read(coadd_fn_new, ext='MASK')
    ivar_new = fitsio.read(coadd_fn_new, ext='IVAR')
    flux_old = fitsio.read(coadd_fn_old, ext='FLUX')
    msk_old = fitsio.read(coadd_fn_old, ext='MASK')
    ivar_old = fitsio.read(coadd_fn_old, ext='IVAR')

    bad = (msk_new != 0) | (ivar_new == 0) | (msk_old != 0) | (ivar_old == 0)
    if np.all(bad):
        return None
    flux_new[bad] = np.nan
    flux_old[bad] = np.nan

    for fiber in idx1:
        fig, ax = plt.subplots(figsize=(12, 4))
        ax.plot(wave, flux_old[fiber], alpha=0.5, lw=3, label='old')
        ax.plot(wave, flux_new[fiber], alpha=0.5, lw=3, label='new')
        ax.grid(alpha=0.5)
        ax.set_xlabel('wavelength')
        ax.set_ylabel('flux')
        ax.set_title(coadd_fn_old + '[{}-{}]'.format(camera+str(petal), fiber))
        ax.legend()
        fig.savefig(f'/global/cfs/cdirs/desicollab/users/rongpu/data/spectro/skypca-y5/plots/m4_cframes/{night}-{expid}-{camera}{petal}-{fiber}_flux.png')
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(12, 4))
        ax.plot(wave, flux_old[fiber] * np.sqrt(ivar_old[fiber]), alpha=0.5, lw=3, label='old')
        ax.plot(wave, flux_new[fiber] * np.sqrt(ivar_new[fiber]), alpha=0.5, lw=3, label='new')
        ax.grid(alpha=0.5)
        ax.set_xlabel('wavelength')
        ax.set_ylabel('flux / flux_error')
        ax.set_title(coadd_fn_old + '[{}-{}]'.format(camera+str(petal), fiber))
        ax.legend()
        fig.savefig(f'/global/cfs/cdirs/desicollab/users/rongpu/data/spectro/skypca-y5/plots/m4_cframes/{night}-{expid}-{camera}{petal}-{fiber}_sn.png')
        plt.close(fig)

tasks = [(idx[ii], petal) for ii in range(len(idx)) for petal in range(10)]

with Pool(processes=256) as pool:
    pool.map(process_exposure_petal, tasks)
