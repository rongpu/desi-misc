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


# transformations are only valid for 0.5<BP-RP<3.0 and 13.5<PHOT_G_MEAN_MAG<17.65
coeffs_bprp = {'g': np.array([-0.0269569852, -0.0569120152, 0.2388830858, -0.2504472232,
        0.1000919830, -0.0137063995]),
 'r': np.array([-0.1041775478, 0.1725203068, -0.1785268995, 0.0978524249,
        -0.0257074158, 0.0026038758]),
 'z': np.array([-0.0346933448, 0.0201737358, 0.0054241207, -0.0188843628,
        0.0096540027, -0.0014980102])}

coeffs_g = {'g': np.array([44.7394180505, -12.9950765935, 1.5190853966, -0.0893304621,
        0.0026441633, -0.0000315555]),
 'r': np.array([667.0475676882, -173.0837779218, 17.5399967510, -0.8596327177,
        0.0200454174, -0.0001722149]),
 'z': np.array([1727.9720400353, -364.5131062097, 25.7499690378, -0.4943381013,
        -0.0158129562, 0.0005468214])}

gaia_columns = ['PHOT_G_MEAN_MAG', 'PHOT_BP_MEAN_MAG', 'PHOT_RP_MEAN_MAG', 'flux_g', 'flux_r', 'flux_z']
sweep_columns = ['RA', 'DEC', 'MASKBITS', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'ANYMASK_G', 'ANYMASK_R', 'ANYMASK_Z']

fns = sorted(glob.glob('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/dr9_north_cross_match/*-gaia.fits'))


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

    mask = (cat['FLUX_G']>0) & (cat['FLUX_R']>0) & (cat['FLUX_Z']>0)
    mask &= (cat['FLUX_IVAR_G']>0) & (cat['FLUX_IVAR_R']>0) & (cat['FLUX_IVAR_Z']>0)
    mask &= (cat['FRACFLUX_G']<0.1) & (cat['FRACFLUX_R']<0.1) & (cat['FRACFLUX_Z']<0.1)
    gaia = gaia[mask]
    cat = cat[mask]
    # print(len(gaia))

    return gaia, cat


n_process = 64
with Pool(processes=n_process) as pool:
    res = pool.map(read_file, fns, chunksize=1)

gaia = vstack([res[ii][0] for ii in range(len(res))])
cat = vstack([res[ii][1] for ii in range(len(res))])
print(len(cat))

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    cat['gmag'] = 22.5 - 2.5*np.log10(cat['FLUX_G'])
    cat['rmag'] = 22.5 - 2.5*np.log10(cat['FLUX_R'])
    cat['zmag'] = 22.5 - 2.5*np.log10(cat['FLUX_Z'])
    gaia['gmag'] = 22.5 - 2.5*np.log10(gaia['flux_g'])
    gaia['rmag'] = 22.5 - 2.5*np.log10(gaia['flux_r'])
    gaia['zmag'] = 22.5 - 2.5*np.log10(gaia['flux_z'])

for band in ['g', 'r', 'z']:
    bprp = (gaia['PHOT_BP_MEAN_MAG']-gaia['PHOT_RP_MEAN_MAG'])
    gaia_g = gaia['PHOT_G_MEAN_MAG'].copy()
    gaia[band+'mag_std'] = gaia[band+'mag'] + poly_val1d(bprp, coeffs_bprp[band])
    gaia[band+'mag_std'] += poly_val1d(gaia_g, coeffs_g[band])

new = Table()
new['RA'] = cat['RA']
new['DEC'] = cat['DEC']
mask0 = (gaia['PHOT_BP_MEAN_MAG']-gaia['PHOT_RP_MEAN_MAG'])>0.6
mask0 &= (gaia['PHOT_BP_MEAN_MAG']-gaia['PHOT_RP_MEAN_MAG'])<2.5
mask0 &= (gaia['PHOT_G_MEAN_MAG']>13.5)
mask0 &= (gaia['PHOT_G_MEAN_MAG']<17.65)

for band in ['g', 'r', 'z']:
    new[band+'mag_diff'] = cat[band+'mag'] - gaia[band+'mag_std']
    new[band+'mag_diff'] = cat[band+'mag'] - gaia[band+'mag_std']
    new[band+'_valid'] = (cat['ANYMASK_'+band.upper()]==0)
    new[band+'_valid'] &= mask0

new.write('/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/misc/gaia_dr9_north_offsets.fits')
