from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

import healpy as hp
# from astropy import units as u
# from astropy.coordinates import SkyCoord

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


gaia_dir = '/project/projectdirs/cosmo/data/gaia/dr2/healpix'
nside = 32
npix = hp.nside2npix(nside)
print('Healpix resolution (arcmin):', np.sqrt(hp.nside2resol(nside, arcmin=True)))

sdss = fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/parent/sdss-specObj-dr16-unique-trimmed.fits')
sdss = Table(sdss)
sdss_hp_idx = hp.ang2pix(nside, sdss['PLUG_RA'], sdss['PLUG_DEC'], nest=True, lonlat=True)
sdss_hp_idx_unique = np.unique(sdss_hp_idx)
print(len(sdss_hp_idx_unique))

sdss_stack = []
gaia_stack = []

for index, hp_idx in enumerate(sdss_hp_idx_unique):

    print(index)

    mask = sdss_hp_idx==hp_idx
    
    gaia_fn = (5-len(str(hp_idx)))*'0'+str(hp_idx)
    gaia = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)))
    gaia = Table(gaia)

    idx1, idx2, d2d, d_ra, d_dec = match_coord(sdss['PLUG_RA'][mask], sdss['PLUG_DEC'][mask], gaia['RA'], gaia['DEC'], search_radius=1., plot_q=False)

    if len(idx1)>0:
        tmp = (sdss[mask])[idx1].copy()
        sdss_stack.append(tmp)
        tmp = gaia[idx2].copy()
        gaia_stack.append(tmp)

print('Done matching!!!!!')

sdss_stack = vstack(sdss_stack)
gaia_stack = vstack(gaia_stack)

sdss_stack.write('/global/cscratch1/sd/rongpu/misc/sdss_gaia_match/sdss.fits')
gaia_stack.write('/global/cscratch1/sd/rongpu/misc/sdss_gaia_match/gaia.fits')

