from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

for fn in ['/pscratch/sd/r/rongpu/ibis/tractor_xmm_stars.fits', '/pscratch/sd/r/rongpu/ibis/tractor_cosmos_stars.fits']:

    # cat = Table(fitsio.read('/pscratch/sd/r/rongpu/ibis/tractor_xmm_stars.fits'))
    cat = Table(fitsio.read(fn))

    nside = 32
    pix = np.unique(hp.ang2pix(nside, cat['ra'], cat['dec'], nest=True, lonlat=True))

    gaia_dir = '/global/cfs/cdirs/cosmo/data/gaia/dr3/healpix'

    gaia = []
    for hp_idx in pix:
        gaia_fn = (5-len(str(hp_idx)))*'0'+str(hp_idx)
        tmp = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), columns=['PHOT_G_MEAN_MAG'])
        # select GAIA_G<13 objects
        idx = np.where(tmp['PHOT_G_MEAN_MAG']<13)[0]
        tmp = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), rows=idx)
        tmp = Table(tmp)
        gaia.append(tmp)
    gaia = vstack(gaia)
    print(len(gaia))

    # convert column names to lower case
    for col in gaia.colnames:
        gaia.rename_column(col, col.lower())

    gaia['bp_rp'] = gaia['phot_bp_mean_mag'] - gaia['phot_rp_mean_mag']
    gaia['pm'] = np.sqrt(gaia['pmra']**2 + gaia['pmdec']**2)

    # Remove duplicates
    if len(np.unique(gaia['source_id']))<len(gaia):
        print('Duplicates exist!')
        gaia.sort('source_id')
        mask = gaia['source_id'][1:]==gaia['source_id'][:-1]
        mask = np.concatenate([[False], mask])
        gaia = gaia[~mask]
        
    # Remove duplicated_source==True
    mask = gaia['duplicated_source'].copy()
    if np.sum(mask)>0:
        print('{} objects with duplicated_source==True'.format(np.sum(mask)))
        gaia = gaia[~mask]
        print(len(gaia))

    # Remove objects with invalid GAIA color
    mask = (~np.isfinite(gaia['bp_rp'])) | (~np.isfinite(gaia['phot_g_mean_mag']))
    mask |= (gaia['phot_g_mean_mag']==0) | (gaia['phot_bp_mean_mag']==0) | (gaia['phot_rp_mean_mag']==0)
    if np.sum(mask)>0:
        print('{} objects with invalid bp_rp'.format(np.sum(mask)))
        gaia = gaia[~mask]
        print(len(gaia))

    print('astrometric error cut:')
    gaia['pm'] = np.sqrt(gaia['pmra']**2 + gaia['pmdec']**2)
    gaia['pmerr'] = np.sqrt(gaia['pmra_error']**2 + gaia['pmdec_error']**2)
    mask = gaia['pm']<30
    print(np.sum(mask)/len(mask))
    mask &= gaia['pmerr']<15
    print(np.sum(mask)/len(mask))
    # mask &= gaia['astrometric_excess_noise']==0
    # print(np.sum(mask)/len(mask))
    # plt.hist(gaia['phot_g_mean_mag'], 100, range=(5, 14), log=True, alpha=0.5)
    # plt.hist(gaia['phot_g_mean_mag'][mask], 100, range=(5, 14), log=True, alpha=0.5)
    # plt.show()

    gaia = gaia[mask]
    print(len(gaia))

    # plt.figure(figsize=(10, 5))
    # plt.plot(gaia['ra'][::2], gaia['dec'][::2], '.', ms=0.5)
    # plt.show()

    ######################## Load GAIA catalog ########################

    # Apply GAIA-LS transformation

    from scipy.interpolate import interp1d

    coeffs = {
    'm411': np.array([-1.21135807e-01,  2.14202972e-01, -4.32275297e-01,  6.33517298e+00,
            -6.16359698e+00, -1.05837512e+01,  3.06209588e+01, -3.23682815e+01,
             1.89619890e+01, -6.68321801e+00,  1.41269689e+00, -1.65375832e-01,
             8.25599733e-03]),
     'm438': np.array([-3.66575264e-02,  1.48435846e-01, -1.73012913e+00,  7.94291264e+00,
            -2.59147985e+00, -1.99128805e+01,  3.67698794e+01, -3.15884298e+01,
             1.59377805e+01, -4.98555542e+00,  9.54462291e-01, -1.02772045e-01,
             4.77913090e-03]),
     'm464': np.array([-1.59291847e-01,  4.04280833e-01,  2.16697525e-01,  1.33091203e+00,
            -1.34037650e+00, -3.65518804e+00,  8.69839989e+00, -8.05791299e+00,
             4.14722148e+00, -1.28622435e+00,  2.40118347e-01, -2.49686426e-02,
             1.11577200e-03]),
     'm490': np.array([-2.88551549e-02,  1.02359685e-01,  4.13203971e-01,  1.58167708e+00,
            -3.59089787e+00, -8.15820939e-01,  9.10026694e+00, -1.13422226e+01,
             7.04242478e+00, -2.53797638e+00,  5.39870363e-01, -6.31013162e-02,
             3.13288100e-03]),
     'm517': np.array([-1.18259758e-01,  3.28822667e-01,  1.34465039e+00, -2.35947092e+00,
            -1.99802997e+00,  7.72242161e+00, -5.60664076e+00, -2.65119231e-01,
             2.37946742e+00, -1.39398479e+00,  3.83994131e-01, -5.32835762e-02,
             3.00229245e-03])
    }

    bprp_min, bprp_max = -0.3, 3.7
    gaia['bp_rp_clipped'] = np.clip(gaia['bp_rp'], bprp_min, bprp_max)

    bands = ['m411', 'm438', 'm464', 'm490', 'm517'] 
    for i, b in enumerate(bands):
        mag = np.copy(gaia['phot_g_mean_mag'])
        for order, c in enumerate(coeffs[b]):
            mag += c * (gaia['bp_rp_clipped'])**order
        gaia[b] = mag

    # gaia.write('/pscratch/sd/r/rongpu/ibis/tractor_xmm_stars-gaia.fits', overwrite=True)
    gaia.write(fn.replace('.fits', '-gaia.fits'), overwrite=False)
