from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
from scipy import stats


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
# mask = tiles['TILERA']<32
# print(np.sum(mask), np.sum(mask)/len(mask))
# tiles = tiles[mask]

#######################################################################


def get_data(tileid, plot=False):
    index = np.where(tiles['TILEID']==tileid)[0][0]
    lastnight, night, expid = tiles['LASTNIGHT'][index], tiles['NIGHT'][index], tiles['EXPID'][index]
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

    if plot:
        plt.figure(figsize=(7, 7))
        plt.plot(cat['TARGET_RA'][~mask], cat['TARGET_DEC'][~mask], '.', ms=2, alpha=1)
        plt.plot(cat['TARGET_RA'][mask], cat['TARGET_DEC'][mask], 'r.', ms=2, alpha=1)
        plt.grid(alpha=0.5)
        # plt.axis([ramax, ramin, decmin, decmax])
        plt.gca().invert_xaxis()
        plt.show()

    cat = cat[mask]
    
    return cat


tileid = 1494
print('TILE:', tileid)
cat = get_data(tileid)
print(len(cat))

mask = cat['psf_corr']<10
print(np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]



