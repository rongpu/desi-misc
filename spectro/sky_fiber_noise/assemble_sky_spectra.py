# srun -N 1 -C cpu -c 256 -t 04:00:00 -q interactive python assemble_sky_spectra.py

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
from multiprocessing import Pool
import itertools


stack_type = 'same_night'
# stack_type = 'different nights 3 months'


def find_sky_fibers(tileid):

    index = np.where(tiles['TILEID']==tileid)[0][0]
    lastnight, night, expid = tiles['LASTNIGHT'][index], tiles['NIGHT'][index], tiles['EXPID'][index]
    dir_path = '/dvs_ro/cfs/cdirs/desi/spectro/redux/loa/tiles/cumulative/{}/*/'.format(tileid)
    fns = glob.glob(os.path.join(dir_path, 'redrock*.fits'))
    fns.sort()
    # print(len(fns), fns[0])

    cat = []
    for fn in fns:
        tmp1 = Table(fitsio.read(fn, ext='FIBERMAP'))
        petal = tmp1['PETAL_LOC'][0]
        expid_str = str(expid).zfill(8)
        # FIBERCORR of each object is identical in the three cameras
        fluxcalib_fn = f'/dvs_ro/cfs/cdirs/desi/spectro/redux/loa/exposures/{night}/{expid_str}/fluxcalib-r{petal}-{expid_str}.fits.gz'
        tmp2 = Table(fitsio.read(fluxcalib_fn, ext='FIBERCORR'))

        mask = tmp1['OBJTYPE']!='SKY'
        median_corr = np.median(tmp2['FLAT_TO_PSF_FLUX'][mask])
        tmp2['psf_corr'] = tmp2['FLAT_TO_PSF_FLUX']/median_corr

        # tmp3 = Table(fitsio.read(fn, ext='REDSHIFTS'))
        # tmp3.remove_columns(['TARGETID'])
        cat.append(hstack([tmp1, tmp2], join_type='exact'))
    cat = vstack(cat)
    cat['LASTNIGHT'] = lastnight
    # cat['LASTNIGHT'] = np.array(os.path.basename(os.path.dirname(fn)), dtype=int)

    mask = cat['OBJTYPE']=='SKY'
    # STUCKPOSITIONER, 1, "INFO: Stuck positioner (but could still be on a good sky location)"
    # POORPOSITION,   10, "Fiber >30 microns from target location"
    mask &= (cat['COADD_FIBERSTATUS']==0) | (cat['COADD_FIBERSTATUS']==2**1) | (cat['COADD_FIBERSTATUS']==(2**1+2**10))
    mask &= cat['psf_corr']<10
    cat = cat[mask]
    
    return cat


tiles = Table.read('/dvs_ro/cfs/cdirs/desi/spectro/redux/loa/tiles-loa.csv')
print(len(tiles))
# tiles['EFFTIME_SPEC_hr'] = tiles['EFFTIME_SPEC'] / 3600.
# tiles['EXPTIME_hr'] = tiles['EXPTIME'] / 3600.
# tiles.sort('EFFTIME_SPEC', reverse=True)
# tiles[:20]

mask = tiles['FAFLAVOR']=='maindark'
tiles = tiles[mask]
print(len(tiles))

mask = tiles['NEXP']==1
tiles = tiles[mask]
print(len(tiles))

# Add per-exposure info
exp = Table.read('/dvs_ro/cfs/cdirs/desi/spectro/redux/loa/exposures-loa.csv')
columns_to_remove = [col for col in exp.colnames if col in tiles.colnames]
columns_to_remove.remove('TILEID')
exp.remove_columns(columns_to_remove)

print(len(tiles))
tiles = join(tiles, exp, keys='TILEID')
print(len(tiles))
tiles.sort('EXPID')

# plt.plot(tiles['SKY_MAG_R_SPEC']-tiles['SKY_MAG_Z_SPEC'], tiles['SKY_MAG_G_SPEC']-tiles['SKY_MAG_R_SPEC'], '.', ms=1)
# plt.xlabel('SKY_MAG_R_SPEC - SKY_MAG_Z_SPEC')
# plt.ylabel('SKY_MAG_G_SPEC - SKY_MAG_R_SPEC')
# plt.show()

# plt.plot(tiles['SKY_MAG_G_SPEC']-tiles['SKY_MAG_R_SPEC'], tiles['SKY_MAG_G_SPEC'], '.', ms=1)
# plt.xlabel('SKY_MAG_G_SPEC - SKY_MAG_R_SPEC')
# plt.ylabel('SKY_MAG_G_SPEC')
# plt.axvline(0.5, lw=1, ls='--', color='r')
# plt.axhline(21.5, lw=1, ls='--', color='r')
# plt.show()

# Pick tiles in moonless conditions
mask = tiles['SKY_MAG_G_SPEC']-tiles['SKY_MAG_R_SPEC'] > 0.5
print(np.sum(mask), np.sum(mask)/len(mask))
mask &= tiles['SKY_MAG_G_SPEC']>21.5
print(np.sum(mask), np.sum(mask)/len(mask))
tiles = tiles[mask]

# # Pick exposures in the Southern imaging (better for identifying blank sky locations)
# mask = tiles['TILEDEC']<30
# print(np.sum(mask), np.sum(mask)/len(mask))
# tiles = tiles[mask]

####################################################################################################
np.random.seed(999)

if stack_type == 'same_night':
    night = 20240212
    mask = tiles['NIGHT']==night
    tiles = tiles[mask]
elif stack_type == 'different nights 3 months':
    mask = (tiles['NIGHT']>20240000) & (tiles['NIGHT']<20240400)
    tiles = tiles[mask]
    tiles = tiles[np.random.choice(len(tiles), size=len(tiles), replace=False)]  # shuffle
    _, idx = np.unique(tiles['NIGHT'], return_index=True)
    tiles = tiles[idx]
    if len(tiles)>50:
        idx = np.random.choice(len(tiles), size=50, replace=False)
        tiles = tiles[idx]
else:
    raise ValueError
####################################################################################################
print('{} tiles:'.format(len(tiles)), list(np.sort(tiles['TILEID'])))

# cat_stack = []
# for index in idx:
#     tileid = tiles['TILEID'][index]
#     cat = find_sky_fibers(tileid)
#     print('TILE {}: {}'.format(tileid, len(cat)))
#     cat_stack.append(cat)
# cat = vstack(cat_stack)
# print(len(cat))

n_process = 32
with Pool(processes=n_process) as pool:
    res = pool.map(find_sky_fibers, tiles['TILEID'])
cat = vstack(res)

fibercount = Table()
fibercount['FIBER'], fibercount['count'] = np.unique(cat['FIBER'], return_counts=True)
fibercount.sort('count', reverse=True)

# select the 500 fibers with the most sky obesrvations
fibercount = fibercount[:500]
fibercount.sort('FIBER')
fibers = list(fibercount['FIBER'])
print(fibers)

mask = np.in1d(cat['FIBER'], fibers)
cat = cat[mask]
print(len(mask), len(cat))

# flux = {'B': {}, 'R': {}, 'Z': {}}
# ivar = {'B': {}, 'R': {}, 'Z': {}}
# msk = {'B': {}, 'R': {}, 'Z': {}}
# resolution = {'B': {}, 'R': {}, 'Z': {}}
# exptime = {'B': {}, 'R': {}, 'Z': {}}

# for camera in ['B', 'R', 'Z']:
#     for fiber in fibers:
#         flux[camera][fiber] = []
#         ivar[camera][fiber] = []
#         msk[camera][fiber] = []
#         resolution[camera][fiber] = []
#         exptime[camera][fiber] = []


def get_spectra(tileid, petal):
    # print(tileid, index, '/', len(tileids))

    mask = (cat['TILEID']==tileid) & (cat['PETAL_LOC']==petal)
    if np.sum(mask)==0:
        return None

    lastnight = cat['LASTNIGHT'][mask][0]
    petal_fibers = cat['FIBER'][mask].copy()

    coadd_fn = f'/dvs_ro/cfs/cdirs/desi/spectro/redux/loa/tiles/cumulative/{tileid}/{lastnight}/coadd-{petal}-{tileid}-thru{lastnight}.fits'
    # print(coadd_fn)

    # fiber, psf_corr = cat['FIBER'][mask1][index], cat['psf_corr'][mask1][index]
    # coadd_index = fiber%500

    fibermap = Table(fitsio.read(coadd_fn, ext='FIBERMAP'))
    redshifts = Table(fitsio.read(coadd_fn.replace('/coadd-', '/redrock-'), ext='REDSHIFTS'))
    assert np.all(fibermap['TARGETID']==fibermap['TARGETID'])
    redshifts = redshifts[['Z', 'ZERR', 'ZWARN', 'CHI2', 'SPECTYPE', 'DELTACHI2']]
    fibermap = hstack([fibermap, redshifts], join_type='exact')
    idx = np.where(np.in1d(fibermap['FIBER'], petal_fibers))[0]
    fibermap = fibermap[idx]

    spec_all = []

    ff_all, ii_all, mm_all, rr_all = {}, {}, {}, {}
    for camera in ['B', 'R', 'Z']:
        ff_all[camera] = fitsio.read(coadd_fn, ext=camera+'_FLUX')[idx]
        ii_all[camera] = fitsio.read(coadd_fn, ext=camera+'_IVAR')[idx]
        mm_all[camera] = fitsio.read(coadd_fn, ext=camera+'_MASK')[idx]==0
        rr_all[camera] = fitsio.read(coadd_fn, ext=camera+'_RESOLUTION')[idx]

    for index, fiber in enumerate(fibermap['FIBER']):

        cat_index = np.where((cat['TILEID']==tileid) & (cat['FIBER']==fiber))[0][0]

        spec = Table()
        spec['FIBER'] = [fiber]
        spec['TILEID'] = [tileid]
        spec['COADD_EXPTIME'] = cat['COADD_EXPTIME'][cat_index]
        spec['COADD_NUMEXP'] = cat['COADD_NUMEXP'][cat_index]
        spec['LASTNIGHT'] = cat['LASTNIGHT'][cat_index]
        psf_corr = cat['psf_corr'][cat_index]
        spec['psf_corr'] = [psf_corr]
        spec['Z'] = [fibermap['Z'][index]]
        spec['ZERR'] = [fibermap['ZERR'][index]]
        spec['SPECTYPE'] = [fibermap['SPECTYPE'][index]]
        spec['DELTACHI2'] = [fibermap['DELTACHI2'][index]]

        for camera in ['B', 'R', 'Z']:
            ff, ii, mm, rr = ff_all[camera][index], ii_all[camera][index], mm_all[camera][index], rr_all[camera][index]
            mm &= ii > 1e-10
            if np.sum(mm)==0:
                print('tile {} fiber {}: all pixels are masked'.format(tileid, fiber))
                continue
            if np.sum(ii[mm] < 1e-10)/np.sum(mm)>0.01:
                print('tile {}: more than 1\% pixels with low IVAR'.format(tileid))
                continue

            ff /= psf_corr
            ii *= psf_corr**2
            ff[~mm] = 0.
            ii[~mm] = 0.

            spec[camera+'_FLUX'] = [ff]
            spec[camera+'_IVAR'] = [ii]
            spec[camera+'_MSK'] = [mm]
            spec[camera+'_RESOLUTION'] = [rr]

        spec_all.append(spec)

    spec_all = vstack(spec_all).filled(0)

    return spec_all


tileids = np.unique(cat['TILEID'])
print(len(tileids))

n_process = 256
arg_list = list(itertools.product(tileids, np.arange(10)))
with Pool(processes=n_process) as pool:
    res = pool.starmap(get_spectra, arg_list)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

spec = vstack(res).filled(0)
spec.sort('FIBER')

if stack_type == 'same_night':
    spec.write('/pscratch/sd/r/rongpu/tmp/sky_spectra/sky_spectra_same_night_20240212.fits', overwrite=True)
elif stack_type == 'different nights 3 months':
    spec.write('/pscratch/sd/r/rongpu/tmp/sky_spectra/sky_spectra_different_nights_202401_202403.fits', overwrite=True)

