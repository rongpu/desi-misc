from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio
import healpy as hp
from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

# gaia_dir = '/dvs_ro/cfs/cdirs/cosmo/data/gaia/dr2/healpix'
gaia_dir = '/dvs_ro/cfs/cdirs/cosmo/data/gaia/dr3/healpix'
nside = 32
npix = hp.nside2npix(nside)
print('Healpix resolution (arcmin):', np.sqrt(hp.nside2resol(nside, arcmin=True)))

sdss = fitsio.read('/dvs_ro/cfs/cdirs/desi/target/analysis/truth/parent/sdss-specObj-dr16-unique-trimmed.fits')
sdss = Table(sdss)
sdss_hp_idx = hp.ang2pix(nside, sdss['PLUG_RA'], sdss['PLUG_DEC'], nest=True, lonlat=True)
sdss_hp_idx_unique = np.unique(sdss_hp_idx)
print(len(sdss_hp_idx_unique))


def process_healpix(hp_idx):
    mask = sdss_hp_idx == hp_idx
    gaia_fn = str(hp_idx).zfill(5)
    gaia = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)))
    gaia = Table(gaia)
    idx1, idx2, d2d, d_ra, d_dec = match_coord(sdss['PLUG_RA'][mask], sdss['PLUG_DEC'][mask], gaia['RA'], gaia['DEC'], search_radius=1., plot_q=False)
    if len(idx1) > 0:
        return ((sdss[mask])[idx1].copy(), gaia[idx2].copy())
    return None


n_workers = 128
print(f'Launching pool with {n_workers} workers over {len(sdss_hp_idx_unique)} healpix pixels')

with Pool(processes=n_workers) as pool:
    results = pool.map(process_healpix, sdss_hp_idx_unique)

sdss_stack = []
gaia_stack = []
for r in results:
    if r is not None:
        sdss_stack.append(r[0])
        gaia_stack.append(r[1])

print('Done matching!!!!!')
sdss_stack = vstack(sdss_stack)
gaia_stack = vstack(gaia_stack)

# sdss_stack.write('/pscratch/sd/r/rongpu/tmp/sdss_gaia_match-gaia_dr2/sdss.fits')
# gaia_stack.write('/pscratch/sd/r/rongpu/tmp/sdss_gaia_match-gaia_dr2/gaia.fits')
sdss_stack.write('/pscratch/sd/r/rongpu/tmp/sdss_gaia_match-gaia_dr3/sdss.fits')
gaia_stack.write('/pscratch/sd/r/rongpu/tmp/sdss_gaia_match-gaia_dr3/gaia.fits')