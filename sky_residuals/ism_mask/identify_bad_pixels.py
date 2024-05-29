from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

import healpy as hp


nside = 256
npix = hp.nside2npix(nside)

sky_south = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_{}_south.fits'.format(nside)))
sky_south['PHOTSYS'] = 'S'
sky_north = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/sky_resid_map_{}_north.fits'.format(nside)))
sky_north['PHOTSYS'] = 'N'

mask = (sky_north['DEC']>32.375)
sky_north = sky_north[mask]
mask = ~np.in1d(sky_south['HPXPIXEL'], sky_north['HPXPIXEL'])
sky = vstack([sky_north, sky_south[mask]])

nsource_min = 250 * (512/nside)**2
mask = sky['nsource_g']>nsource_min
mask &= sky['nsource_r']>nsource_min
mask &= sky['nsource_z']>nsource_min
print(np.sum(mask)/len(mask))
sky = sky[mask]

mask_bad = (sky['sky_median_g']>0.00065) | (sky['sky_median_g']<-0.0003)
mask_bad |= (sky['sky_median_r']>0.00115) | (sky['sky_median_r']<-0.00045)
mask_bad |= (sky['PHOTSYS']=='S') & ((sky['sky_median_z']<-0.002) | (sky['sky_median_z']>0.0025))
mask_bad |= (sky['PHOTSYS']=='N') & ((sky['sky_median_z']<-0.004) | (sky['sky_median_z']>-0.0007))
print(np.sum(mask_bad)/len(mask_bad))

mask_bad |= (sky['PHOTSYS']=='S') & ((sky['sky_nmad_g']>0.00315) | (sky['sky_nmad_r']>0.00525))
mask_bad |= (sky['PHOTSYS']=='N') & ((sky['sky_nmad_g']>0.0036) | (sky['sky_nmad_r']>0.0071))
mask_bad |= (sky['sky_nmad_z']>0.012)
print(np.sum(mask_bad)/len(mask_bad))

npix = hp.nside2npix(nside)
tmp = np.zeros(npix)
tmp[sky['HPXPIXEL'][mask_bad]] = 1.
tmp = hp.ud_grade(tmp, nside*2, order_in='RING', order_out='RING')
tmp = np.where(tmp!=0)[0]

bad_pixels = np.unique(np.concatenate([tmp, hp.get_all_neighbours(nside*2, tmp).flatten()]))
tt = Table()
tt['HPXPIXEL'] = bad_pixels
mask = np.in1d(tt['HPXPIXEL'], tmp)
tt['CORE'] = mask.copy()
tt.write('/global/cfs/cdirs/desi/users/rongpu/imaging_mc/ism_mask/bad_pixels_v1_{}_ring.fits'.format(nside*2))
