from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from scipy import stats
from multiprocessing import Pool


def get_data_with_fibercorr_factor(index):
    tileid, lastnight, night, expid = tiles['TILEID'][index], tiles['LASTNIGHT'][index], tiles['NIGHT'][index], tiles['EXPID'][index]
    dir_path = '/global/cfs/cdirs/desi/spectro/redux/loa/tiles/cumulative/{}/*/'.format(tileid)
    fns = glob.glob(os.path.join(dir_path, 'redrock*.fits'))
    fns.sort()
    # print(len(fns), fns[0])
    cat = []
    for fn in fns:
        tmp1 = Table(fitsio.read(fn, ext='FIBERMAP'))
        petal = tmp1['PETAL_LOC'][0]
        expid_str = str(expid).zfill(8)
        # FIBERCORR of each object is identical in the three cameras
        fluxcalib_fn = f'/global/cfs/cdirs/desi/spectro/redux/loa/exposures/{night}/{expid_str}/fluxcalib-r{petal}-{expid_str}.fits.gz'
        tmp2 = Table(fitsio.read(fluxcalib_fn, ext='FIBERCORR'))
        cat.append(hstack([tmp1, tmp2], join_type='exact'))
    cat = vstack(cat)
    cat['LASTNIGHT'] = lastnight
    # cat['LASTNIGHT'] = np.array(os.path.basename(os.path.dirname(fn)), dtype=int)
    return cat


tiles = Table.read('/global/cfs/cdirs/desi/spectro/redux/loa/tiles-loa.csv')
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
exp = Table.read('/global/cfs/cdirs/desi/spectro/redux/loa/exposures-loa.csv')
columns_to_remove = [col for col in exp.colnames if col in tiles.colnames]
columns_to_remove.remove('TILEID')
exp.remove_columns(columns_to_remove)

print(len(tiles))
tiles = join(tiles, exp, keys='TILEID')
print(len(tiles))
tiles.sort('EXPID')

plt.plot(tiles['SKY_MAG_R_SPEC']-tiles['SKY_MAG_Z_SPEC'], tiles['SKY_MAG_G_SPEC']-tiles['SKY_MAG_R_SPEC'], '.', ms=1)
plt.xlabel('SKY_MAG_R_SPEC - SKY_MAG_Z_SPEC')
plt.ylabel('SKY_MAG_G_SPEC - SKY_MAG_R_SPEC')
plt.show()

plt.plot(tiles['SKY_MAG_G_SPEC']-tiles['SKY_MAG_R_SPEC'], tiles['SKY_MAG_G_SPEC'], '.', ms=1)
plt.xlabel('SKY_MAG_G_SPEC - SKY_MAG_R_SPEC')
plt.ylabel('SKY_MAG_G_SPEC')
plt.axvline(0.5, lw=1, ls='--', color='r')
plt.axhline(21.5, lw=1, ls='--', color='r')
plt.show()

# Pick tiles in moonless conditions
mask = tiles['SKY_MAG_G_SPEC']-tiles['SKY_MAG_R_SPEC'] > 0.5
print(np.sum(mask), np.sum(mask)/len(mask))
mask &= tiles['SKY_MAG_G_SPEC']>21.5
print(np.sum(mask), np.sum(mask)/len(mask))
tiles = tiles[mask]

np.random.seed(888)
idx = np.random.choice(len(tiles), size=100, replace=False)
n_process = 128
with Pool(processes=n_process) as pool:
    res = pool.map(get_data_with_fibercorr_factor, idx)
cat = vstack(res, join_type='exact')

cat.write('/pscratch/sd/r/rongpu/tmp/fiberapcorr.fits')


