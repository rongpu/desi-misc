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
# field = 'north'

# sweep_dir = '/global/cscratch1/sd/adamyers/dr9m-sep26-2020/{}/sweep'.format(field)
sweep_dir = '/global/cscratch1/sd/landriau/dr9m-partial-sweep/{}/sweep'.format(field)
tractor_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr9m/{}/tractor'.format(field)
output_dir = '/global/cscratch1/sd/rongpu/desi/dr9m-oct26-2020_sweep_apflux_resid/{}'.format(field)

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
    data_init = list(-99.*np.ones((6, len(sweep), 8), dtype='float32'))
    apflux_resid = Table(data=data_init, names=['apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z', 'apflux_blobresid_g', 'apflux_blobresid_r', 'apflux_blobresid_z'])
    for colname in ['apflux_resid_w1', 'apflux_resid_w2']:
        apflux_resid[colname] = -99 * np.ones((len(sweep), 5), dtype='float32')

    for brickname in all_brickname:
        # print(brickname)
        brick_mask = sweep['BRICKNAME']==brickname

        tractor_path = os.path.join(tractor_dir, brickname[:3], 'tractor-'+brickname+'.fits')
        tractor = fitsio.read(tractor_path, columns=['objid', 'apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z', 'apflux_blobresid_g', 'apflux_blobresid_r', 'apflux_blobresid_z', 'apflux_resid_w1', 'apflux_resid_w2'])
        tractor = Table(tractor)  # this resolves the problem with numpy.bytes_
        
        ########## SANITY CHECKS ##########
        # Two assumptions (for both sweep and tractor catalogs):
        # 1. OBJID values for each brick are unique
        # 2. OBJID values for each brick are monotonically increasing
        if not (len(np.unique(tractor['objid']))==len(tractor)):
            raise ValueError('tractor OBJID values are not unique! sweep_index = {}'.format(sweep_index))
        if not (np.all(np.diff(tractor['objid'])>0)):
            raise ValueError('tractor OBJID values are monotonically increasing! sweep_index = {}'.format(sweep_index))
        if not (len(np.unique(sweep['OBJID'][brick_mask]))==np.sum(brick_mask)):
            raise ValueError('sweep OBJID values are not unique! sweep_index = {}'.format(sweep_index))
        if not (np.all(np.diff(sweep['OBJID'][brick_mask])>0)):
            raise ValueError('sweep OBJID values are monotonically increasing! sweep_index = {}'.format(sweep_index))
        ###################################
        
        mask = np.in1d(tractor['objid'], sweep['OBJID'][brick_mask])
        # another sanity check
        if not np.all(tractor['objid'][mask]==sweep['OBJID'][brick_mask]):
            ###############################################################################################################
            ##################################### TEMPORARY!!!!! Only for partial sweeps ##################################
            ##################################### Restore ValueError in full release!!!!! #################################
            ###############################################################################################################
            print('OBJID does not match! sweep_index = {}, brickname = {}'.format(sweep_index, brickname))
            continue
            # raise ValueError('OBJID does not match! sweep_index = {}'.format(sweep_index))
            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
            ###############################################################################################################
        
        for colname in ['apflux_resid_g', 'apflux_resid_r', 'apflux_resid_z', 'apflux_blobresid_g', 'apflux_blobresid_r', 'apflux_blobresid_z', 'apflux_resid_w1', 'apflux_resid_w2']:
            apflux_resid[colname][brick_mask] = tractor[colname][mask]

    # clear cache
    gc.collect()

    output_path = os.path.join(output_dir, sweep_fn[:-5]+'-resid.fits')
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
