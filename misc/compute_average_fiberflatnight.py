from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
from multiprocessing import Pool

nights = np.sort([int(tmp.split('/')[-1]) for tmp in glob.glob('/dvs_ro/cfs/cdirs/desi/spectro/redux/daily/calibnight/*')])
# mask = nights > 20230700
# nights = nights[mask]
print(nights)
print(len(nights))

# fiberxy = []
# for petal in range(10):
#     fiberxy.append(Table(fitsio.read(f'/dvs_ro/cfs/cdirs/desi/spectro/redux/daily/calibnight/20231011/fiberflatnight-r{petal}-20231011.fits', ext='FIBERMAP')))
# fiberxy = vstack(fiberxy)[['FIBER', 'FIBERASSIGN_X', 'FIBERASSIGN_Y']]
# assert np.all(fiberxy['FIBER']==np.arange(5000))


def average_fiberflat(night):
    if len(glob.glob(f'/dvs_ro/cfs/cdirs/desi/spectro/redux/daily/calibnight/{night}/fiberflatnight-*'))==0:
        return None
    cat = Table()
    cat['NIGHT'] = night*np.ones(5000, dtype='int32')
    cat['FIBER'] = np.arange(5000, dtype='int32')
    cat['B_MEDIAN'], cat['B_MEAN'], cat['R_MEDIAN'], cat['R_MEAN'], cat['Z_MEDIAN'], cat['Z_MEAN'] = -99*np.ones((6, 5000), dtype='float32')
    cat['FIBERSTATUS'] = -99*np.ones(5000, dtype=int)
    for band in ['b', 'r', 'z']:
        for petal in range(10):
            fn = f'/dvs_ro/cfs/cdirs/desi/spectro/redux/daily/calibnight/{night}/fiberflatnight-{band}{petal}-{night}.fits'
            if not os.path.isfile(fn):
                continue
            ff = fitsio.read(fn, ext='FIBERFLAT')
            assert len(ff)==500
            cat[band.upper()+'_MEDIAN'][petal*500:petal*500+500] = np.median(ff, axis=1)
            cat[band.upper()+'_MEAN'][petal*500:petal*500+500] = np.mean(ff, axis=1)
            if band=='r':
                cat['FIBERSTATUS'][petal*500:petal*500+500] = fitsio.read(fn, ext='FIBERMAP')['FIBERSTATUS']
    return cat


n_process = 256
with Pool(processes=n_process) as pool:
    res = pool.map(average_fiberflat, nights)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res)
print(len(cat))
cat.write('/pscratch/sd/r/rongpu/tmp/fiberflatnight.fits', overwrite=True)


############################################## Get Parallactic and dome DOMEAZ angles ##############################################

def read_header(night):
    fns = glob.glob(f'/dvs_ro/cfs/cdirs/desi/spectro/redux/daily/calibnight/{night}/fiberflatnight-*')
    fns.sort()
    if len(fns)==0:
        return None
    cat = Table()
    cat['NIGHT'] = [night]
    fn = fns[0]
    header = fitsio.read_header(fn, ext=0)

    if 'PARALLAC' in header.keys():
        cat['PARALLAC'] = [header['PARALLAC']]  # [deg] Parallactic angle
    else:
        cat['PARALLAC'] = -99.

    if 'DOMEAZ' in header.keys():
        cat['DOMEAZ'] = [header['DOMEAZ']]  # [deg] Dome azimuth angle
    else:
        cat['DOMEAZ'] = -99.

    if 'MJD-OBS' in header.keys():
        cat['MJD-OBS'] = [header['MJD-OBS']]
    else:
        cat['MJD-OBS'] = -99.

    return cat


n_process = 256
with Pool(processes=n_process) as pool:
    res = pool.map(read_header, nights)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res)
print(len(cat))
cat.write('/pscratch/sd/r/rongpu/tmp/fiberflatnight_headers.fits', overwrite=True)
