# Created line-matched aperture flux catalog for the sweeps

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
import gc
from multiprocessing import Pool

n_processes = 32
field = 'south'

sweep_dir = '/global/project/projectdirs/cosmo/data/legacysurvey/dr8/'+field+'/sweep/8.0'
tractor_dir = '/global/project/projectdirs/cosmo/data/legacysurvey/dr8/'+field+'/tractor'

sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]

def get_apflux_resid(sweep_index):

    sweep_fn = sweep_fn_list[sweep_index]
    sweep = fitsio.read(os.path.join(sweep_dir, sweep_fn), columns=['BRICKNAME', 'OBJID'])
    sweep = Table(sweep)  # this resolves the problem with numpy.bytes_
    print(sweep_fn, len(sweep))

    all_brickname = np.unique(sweep['BRICKNAME'])
    print(len(all_brickname))

    # -99 for objects no found in tractor catalogs (which should not happen)
    data_init = list(-99.*np.ones((3, len(sweep), 8), dtype='float32'))
    apflux_resid = Table(data=data_init, names=['apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z'])

    for brickname in all_brickname:
        # print(brickname)
        brick_mask = sweep['BRICKNAME']==brickname

        tractor_path = os.path.join(tractor_dir, brickname[:3], 'tractor-'+brickname+'.fits')
        tractor = fitsio.read(tractor_path, columns=['objid', 'apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z'])
        tractor = Table(tractor)  # this resolves the problem with numpy.bytes_
        
        ########## SANITY CHECKS ##########
        # Two assumptions (for both sweep and tractor catalogs):
        # 1. OBJID values for each brick are unique
        # 2. OBJID values for each brick are monotonically increasing
        if not (len(np.unique(tractor['objid']))==len(tractor)):
            raise ValueError('tractor OBJID values are not unique!')
        if not (np.all(np.diff(tractor['objid'])>0)):
            raise ValueError('tractor OBJID values are monotonically increasing!')
        if not (len(np.unique(sweep['OBJID'][brick_mask]))==np.sum(brick_mask)):
            raise ValueError('sweep OBJID values are not unique!')
        if not (np.all(np.diff(sweep['OBJID'][brick_mask])>0)):
            raise ValueError('sweep OBJID values are monotonically increasing!')
        ###################################
        
        mask = np.in1d(tractor['objid'], sweep['OBJID'][brick_mask])
        # another sanity check
        if not np.all(tractor['objid'][mask]==sweep['OBJID'][brick_mask]):
            raise ValueError('OBJID does not match!')
        
        for colname in ['apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z']:
            apflux_resid[colname][brick_mask] = tractor[colname][mask]

    # clear cache
    gc.collect()

    output_path = os.path.join('/global/cscratch1/sd/rongpu/desi/dr8_sweep_apflux_resid/', field, sweep_fn[:-5]+'-resid.fits')
    apflux_resid.write(output_path)

    # return apflux_resid
    return None

if __name__ == '__main__':

    print('Start!')

    time_start = time.time()

    # start multiple worker processes
    with Pool(processes=n_processes) as pool:
        pool.map(get_apflux_resid, range(len(sweep_fn_list)))

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
