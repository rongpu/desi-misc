from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

keys = ['SKYLEVEL', 'ETCSKY', 'ETCSKYLV']


def read_header_info(index):
    fn = "/global/cfs/cdirs/desi/spectro/data/{}/{:08d}/desi-{:08d}.fits.fz".format(cat["NIGHT"][index],cat["EXPID"][index],cat["EXPID"][index])
    header = fitsio.read_header(fn, ext=1)
    list = []
    for key in keys:
        list.append(header[key])
    # HUMIDITY, OUTTEMP, DEWPOINT, EWALLCOU
    return np.array(list)


cat = Table(fitsio.read('/global/cfs/cdirs/desi/spectro/redux/daily/exposures-daily.fits'))
mask = cat['SURVEY']=='main'
mask &= cat['PROGRAM']=='dark'
mask &= cat['EFFTIME_ETC']>100
mask &= cat['EFFTIME_SPEC']>100
mask &= cat['MJD'] != 0
mask &= cat['EFFTIME_GFA']!=0
cat = cat[mask]
print(len(cat))

###############################################################
# Earlier exposures do not have all the ETC data in the header
mask = cat['NIGHT']>20210900
cat = cat[mask]
print(len(cat))
###############################################################


from multiprocessing import Pool
n_process = 256
with Pool(processes=n_process) as pool:
    res = pool.map(read_header_info, np.arange(len(cat)), chunksize=1)
res = np.stack(res)

tt = Table()
tt['EXPID'] = cat['EXPID']
for ii, key in enumerate(keys):
    tt[key] = res[:, ii]
tt.write('/pscratch/sd/r/rongpu/tmp/exposures-daily-header-data.fits')
