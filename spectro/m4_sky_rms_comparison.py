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

mask = exposures['NIGHT']>20250700
exposures = exposures[mask]
print(len(exposures))

mask = exposures['EXPTIME']>800
exposures = exposures[mask]
print('EXPTIME', len(exposures))

camera = 'r'

# np.random.seed(999)
# idx = np.sort(np.random.choice(len(exposures), size=20, replace=False))
idx = np.arange(len(exposures))

print(np.unique(exposures['NIGHT'][idx]))

def process_exposure_petal(args):
    index, petal = args
    night = exposures['NIGHT'][index]
    expid = exposures['EXPID'][index]
    expid_str = str(expid).zfill(8)

    coadd_fn_new = f'/dvs_ro/cfs/cdirs/desi/spectro/redux/m4/exposures/{night}/{expid_str}/cframe-{camera}{petal}-{expid_str}.fits.gz'
    if not os.path.isfile(coadd_fn_new):
        print('file does not exist:', coadd_fn_new)
        return

    coadd_fn_old = f'/dvs_ro/cfs/cdirs/desi/spectro/redux/loa/exposures/{night}/{expid_str}/cframe-{camera}{petal}-{expid_str}.fits.gz'
    if not os.path.isfile(coadd_fn_old):
        coadd_fn_old = f'/dvs_ro/cfs/cdirs/desi/spectro/redux/daily/exposures/{night}/{expid_str}/cframe-{camera}{petal}-{expid_str}.fits.gz'

    cat = Table(fitsio.read(coadd_fn_new, ext='FIBERMAP'))
    mask = (cat['OBJTYPE']=='SKY')
    mask &= (cat['FIBERSTATUS']==0) | (cat['FIBERSTATUS']==2**1) | (cat['FIBERSTATUS']==(2**1+2**10))
    idx1 = np.where(mask)[0]

    wave = fitsio.read(coadd_fn_new, ext='WAVELENGTH')
    flux_new = fitsio.read(coadd_fn_new, ext='FLUX')
    msk_new = fitsio.read(coadd_fn_new, ext='MASK')
    ivar_new = fitsio.read(coadd_fn_new, ext='IVAR')
    flux_old = fitsio.read(coadd_fn_old, ext='FLUX')
    msk_old = fitsio.read(coadd_fn_old, ext='MASK')
    ivar_old = fitsio.read(coadd_fn_old, ext='IVAR')

    bad = (msk_new != 0) | (ivar_new == 0) | (msk_old != 0) | (ivar_old == 0)
    flux_new[bad] = np.nan
    flux_old[bad] = np.nan

    rms_ratio = []
    for fiber in idx1:
        if np.sum(~np.isnan(flux_new[fiber]))<100:
            continue
        rms_ratio.append(np.nanstd(flux_new[fiber]) / np.nanstd(flux_old[fiber]))

    if len(rms_ratio)>0:
        tt = Table()
        tt['rms_ratio'] = rms_ratio
        tt['EXPID'] = expid
        tt['NIGHT'] = night
        tt['PETAL'] = petal
        return tt
    else:
        print('No valid sky fiber', night, expid, petal)
        return None

tasks = [(idx[ii], petal) for ii in range(len(idx)) for petal in range(10)]

with Pool(processes=256) as pool:
    res = pool.map(process_exposure_petal, tasks)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

tt = vstack(res)
tt.write('/global/cfs/cdirs/desicollab/users/rongpu/data/spectro/skypca-y5/m4_rms_ratio-{}-more_sky_fibers.fits'.format(camera))
