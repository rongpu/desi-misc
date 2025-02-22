from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool


def poly_val1d(x, m):
    '''
    Evaluate the 1D polynomial from x values and polynomial parameters

    See also rlm_fit1d.py

    INPUT:
    1D array of x values;
    1D array of polynomial parameters (for example generated by rlm_fit1d.py).

    OUTPUT:
    1D array of the evaluated values of the polynomial.
    '''

    order = len(m)-1
    z = np.zeros(x.shape)
    for i in range(order+1):
        z += m[i] * x**i
    return z


coeffs_bprp = {'g': np.array([0.0416867638, -0.2958055870, 0.4922869977, -0.4032130276,
        0.1459382610, -0.0190536024]),
 'r': np.array([-0.0529926024, 0.0361338856, -0.0188637440, -0.0004627745,
        0.0025104619, -0.0004741510]),
 'i': np.array([-0.0873044994, 0.1770173852, -0.1933434841, 0.0996162932,
        -0.0245789844, 0.0023553514]),
 'z': np.array([-0.0431176311, 0.0586222512, -0.0112472116, -0.0188292696,
        0.0109515264, -0.0016498615])}

coeffs_g = {'g': np.array([-0.6417454609, 0.9144340164, -0.1939379255, 0.0165146864,
        -0.0006337951, 0.0000091121]),
 'r': np.array([-17.4056692603, 5.5389908307, -0.7004043218, 0.0440353839,
        -0.0013759207, 0.0000170644]),
 'i': np.array([112.4088645838, -34.6121088481, 4.2534703382, -0.2607416920,
        0.0079739763, -0.0000973534]),
 'z': np.array([50.9842611342, -15.8466786120, 1.9770712887, -0.1236950359,
        0.0038797140, -0.0000488000])}

gaia_columns = ['PHOT_G_MEAN_MAG', 'PHOT_BP_MEAN_MAG', 'PHOT_RP_MEAN_MAG', 'flux_g', 'flux_r', 'flux_i', 'flux_z']
sweep_columns = ['RA', 'DEC', 'MASKBITS', 'FLUX_G', 'FLUX_R', 'FLUX_I', 'FLUX_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_I', 'FLUX_IVAR_Z', 'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_I', 'FRACFLUX_Z', 'ANYMASK_G', 'ANYMASK_R', 'ANYMASK_I', 'ANYMASK_Z']

fns = sorted(glob.glob('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/dr10_south_cross_match/*-gaia.fits'))


def read_file(fn):

    print(fn)
    gaia = Table(fitsio.read(fn, columns=gaia_columns))
    cat = Table(fitsio.read(fn.replace('-gaia.fits', '-ls.fits'), columns=sweep_columns))
    # print(len(gaia))

    mask = gaia['PHOT_G_MEAN_MAG']<17.65
    gaia = gaia[mask]
    cat = cat[mask]
    # print(len(gaia))

    # maskbits = [0, 2, 3, 4, 5, 6, 7, 10, 12, 13]
    maskbits = [0, 1, 10, 12, 13]
    mask_clean = np.ones(len(cat), dtype=bool)
    for bit in maskbits:
        mask_clean &= (cat['MASKBITS'] & 2**bit)==0
    # print(np.sum(~mask_clean)/len(mask_clean))
    gaia = gaia[mask_clean]
    cat = cat[mask_clean]
    # print(len(gaia))

    # mask = (cat['FLUX_G']>0) & (cat['FLUX_R']>0) & (cat['FLUX_Z']>0) & (cat['FLUX_I']>0)
    # mask &= (cat['FLUX_IVAR_G']>0) & (cat['FLUX_IVAR_R']>0) & (cat['FLUX_IVAR_Z']>0) & (cat['FLUX_IVAR_I']>0)
    # mask &= (cat['FRACFLUX_G']<0.1) & (cat['FRACFLUX_R']<0.1) & (cat['FRACFLUX_Z']<0.1) & (cat['FRACFLUX_I']<0.1)
    # gaia = gaia[mask]
    # cat = cat[mask]
    # print(len(gaia))

    return gaia, cat


n_process = 64
with Pool(processes=n_process) as pool:
    res = pool.map(read_file, fns, chunksize=1)

gaia = vstack([res[ii][0] for ii in range(len(res))])
cat = vstack([res[ii][1] for ii in range(len(res))])

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    cat['gmag'] = 22.5 - 2.5*np.log10(cat['FLUX_G'])
    cat['rmag'] = 22.5 - 2.5*np.log10(cat['FLUX_R'])
    cat['imag'] = 22.5 - 2.5*np.log10(cat['FLUX_I'])
    cat['zmag'] = 22.5 - 2.5*np.log10(cat['FLUX_Z'])
    gaia['gmag'] = 22.5 - 2.5*np.log10(gaia['flux_g'])
    gaia['rmag'] = 22.5 - 2.5*np.log10(gaia['flux_r'])
    gaia['imag'] = 22.5 - 2.5*np.log10(gaia['flux_i'])
    gaia['zmag'] = 22.5 - 2.5*np.log10(gaia['flux_z'])

for band in ['g', 'r', 'i', 'z']:
    bprp = (gaia['PHOT_BP_MEAN_MAG']-gaia['PHOT_RP_MEAN_MAG'])
    gaia_g = gaia['PHOT_G_MEAN_MAG'].copy()
    gaia[band+'mag_std'] = gaia[band+'mag'] + poly_val1d(bprp, coeffs_bprp[band])
    gaia[band+'mag_std'] += poly_val1d(gaia_g, coeffs_g[band])

new = Table()
new['RA'] = cat['RA']
new['DEC'] = cat['DEC']
mask0 = (gaia['PHOT_BP_MEAN_MAG']-gaia['PHOT_RP_MEAN_MAG'])>0.6
mask0 &= (gaia['PHOT_BP_MEAN_MAG']-gaia['PHOT_RP_MEAN_MAG'])<2.5

for band in ['g', 'r', 'i', 'z']:
    new[band+'mag_diff'] = cat[band+'mag'] - gaia[band+'mag_std']
    new[band+'mag_diff'] = cat[band+'mag'] - gaia[band+'mag_std']
    new[band+'_valid'] = (cat['ANYMASK_'+band.upper()]==0)
    new[band+'_valid'] &= (cat['FLUX_'+band.upper()]>0) & (cat['FLUX_IVAR_'+band.upper()]>0) & (cat['FRACFLUX_'+band.upper()]<0.1)
    new[band+'_valid'] &= mask0

new.write('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_dr10_offsets.fits')
