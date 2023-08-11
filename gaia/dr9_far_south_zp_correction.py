# Apply ZP correction for objects at DEC<-29.25 and save the corrected sweep catalogs

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

from multiprocessing import Pool


def decode_sweep_name(sweepname):
    # taken from desihub/desitarget

    sweepname = os.path.basename(sweepname)

    ramin, ramax = float(sweepname[6:9]), float(sweepname[14:17])
    decmin, decmax = float(sweepname[10:13]), float(sweepname[18:21])

    if sweepname[9] == 'm':
        decmin *= -1
    if sweepname[17] == 'm':
        decmax *= -1

    return [ramin, ramax, decmin, decmax]


mdpole_dict = {#'g': [-0.006279329661472231, [0.00732357, -0.00037836, -0.00413991]],
               'r': [-0.005438745197063, [0.00253447, -0.00363545, 0.01021869]],
               'z': [0.014450568216766158, [6.41211157e-05, -1.35882444e-02, 3.91261375e-03]]}

sweep_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0'
output_dir = '/pscratch/sd/r/rongpu/dr9_desi_photoz/sweep_zp_corrected'


sweep_all_path = sorted(glob.glob(os.path.join(sweep_dir, '*.fits')))
sweep_fn_list = [os.path.basename(sweep_all_path[ii]) for ii in range(len(sweep_all_path))]
sweep_fn_list_new = []
for sweep_fn in sweep_fn_list:
    if decode_sweep_name(sweep_fn)[3]<=-25:
        sweep_fn_list_new.append(sweep_fn)
sweep_fn_list = sweep_fn_list_new
print(len(sweep_fn_list))


def do_stuff(sweep_fn):

    sweep_path = os.path.join(sweep_dir, sweep_fn)
    cat = Table(fitsio.read(sweep_path))

    # Convert (RA, DEC) to Cartesian coordinates
    x = np.cos(cat['RA']/180*np.pi)*np.cos(cat['DEC']/180*np.pi)
    y = np.sin(cat['RA']/180*np.pi)*np.cos(cat['DEC']/180*np.pi)
    z = np.sin(cat['DEC']/180*np.pi)

    mask_corr = cat['DEC']<-29.25

    for band in ['r', 'z']:

        monopole, dipole = mdpole_dict[band]
        vv = np.array([x, y, z]).T
        md_map = monopole + np.dot(vv, dipole)
        flux_rescale = 10**(md_map*0.4)
        flux_rescale[~mask_corr] = 1.

        cat['FLUX_'+band.upper()] = cat['FLUX_'+band.upper()] * flux_rescale
        cat['FIBERFLUX_'+band.upper()] = cat['FIBERFLUX_'+band.upper()] * flux_rescale
        cat['FIBERTOTFLUX_'+band.upper()] = cat['FIBERTOTFLUX_'+band.upper()] * flux_rescale
        cat['FLUX_IVAR_'+band.upper()] = cat['FLUX_IVAR_'+band.upper()] / flux_rescale**2

    cat.write(os.path.join(output_dir, sweep_fn))

    return None


n_processes = 64
with Pool(processes=n_processes) as pool:
    pool.map(do_stuff, sweep_fn_list)


