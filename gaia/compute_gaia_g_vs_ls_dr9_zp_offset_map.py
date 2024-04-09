from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp
from multiprocessing import Pool

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


gaia_dir = '/global/cfs/cdirs/cosmo/data/gaia/dr3/healpix'
gaia_nside = 32
gaia_npix = hp.nside2npix(gaia_nside)

sweep_fn_list = glob.glob('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/sweep/9.0/*.fits')


def get_sweep_and_gaia(sweep_fn):

    cat = Table(fitsio.read(sweep_fn,
                columns=['RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'TYPE', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FRACFLUX_G', 'FRACFLUX_R', 'FRACFLUX_Z', 'FRACMASKED_G', 'FRACMASKED_R', 'FRACMASKED_Z', 'FIBERFLUX_G', 'FIBERTOTFLUX_G', 'FIBERFLUX_R', 'FIBERTOTFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_Z', 'ANYMASK_G', 'ANYMASK_R', 'ANYMASK_Z', 'REF_EPOCH', 'PMRA', 'PMDEC', 'EBV']))
    mask = 22.5 - 2.5*np.log10(cat['FLUX_R']) < 20.0
    # mask &= cat['TYPE']!='DUP'
    mask &= cat['TYPE']=='PSF'
    mask &= cat['FRACFLUX_R']<0.01
    mask &= (cat['NOBS_G']>0) & (cat['NOBS_R']>0) & (cat['NOBS_Z']>0)
    mask &= (cat['FLUX_IVAR_G']>0) & (cat['FLUX_IVAR_R']>0) & (cat['FLUX_IVAR_Z']>0)
    mask &= (cat['FRACFLUX_G']<0.01) & (cat['FRACFLUX_R']<0.01) & (cat['FRACFLUX_Z']<0.01)
    mask &= (cat['FRACMASKED_G']<0.6) & (cat['FRACMASKED_R']<0.6) & (cat['FRACMASKED_Z']<0.6)
    mask &= (cat['FIBERFLUX_G']/cat['FIBERTOTFLUX_G']>0.99) & (cat['FIBERFLUX_R']/cat['FIBERTOTFLUX_R']>0.99) & (cat['FIBERFLUX_Z']/cat['FIBERTOTFLUX_Z']>0.99)
    mask &= (cat['ANYMASK_G']==0) & (cat['ANYMASK_R']==0) & (cat['ANYMASK_Z']==0)
    cat = cat[mask]

    if len(cat)==0:
        return None

    # The reference epoch for Gaia DR3 (both Gaia EDR3 and the full Gaia DR3) is 2016.0
    # https://www.cosmos.esa.int/web/gaia/dr3
    cat['RA'] = cat['RA'] - (cat['REF_EPOCH'] - 2016.0) * cat['PMRA'] * 1e-3/3600 / np.cos(np.radians(cat['DEC']))
    cat['RA'] = (cat['RA'] + 360)%360  # Wrap around
    cat['DEC'] = cat['DEC'] - (cat['REF_EPOCH'] - 2016.0) * cat['PMDEC'] * 1e-3/3600

    healpix_list = np.unique(hp.ang2pix(gaia_nside, cat['RA'], cat['DEC'], nest=True, lonlat=True))

    gaia = []
    for hp_index in healpix_list:
        gaia_fn = str(hp_index).zfill(5)
        gaia_chunk = Table(fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), columns=['SOURCE_ID', 'RA', 'DEC', 'PHOT_G_MEAN_MAG', 'PHOT_BP_MEAN_MAG', 'PHOT_RP_MEAN_MAG', 'PHOT_G_MEAN_FLUX_OVER_ERROR']))
        mask = (gaia_chunk['PHOT_G_MEAN_MAG']>14) & (gaia_chunk['PHOT_G_MEAN_MAG']<18)
        gaia_chunk = gaia_chunk[mask]
        gaia.append(gaia_chunk)
    gaia = vstack(gaia)
    if len(gaia)==0:
        return None
    idx1, idx2, d2d, d_ra, d_dec = match_coord(cat['RA'], cat['DEC'], gaia['RA'], gaia['DEC'], search_radius=0.1, plot_q=False)

    if len(idx1)==0:
        return None

    cat = cat[idx1]
    gaia = gaia[idx2]
    cat = cat['RA', 'DEC', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'MASKBITS', 'EBV']
    gaia.remove_columns(['RA', 'DEC'])

    cat = hstack([cat, gaia])

    return cat


n_process = 128
with Pool(processes=n_process) as pool:
    res = pool.map(get_sweep_and_gaia, sweep_fn_list, chunksize=1)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

cat = vstack(res)
cat.write('/pscratch/sd/r/rongpu/tmp/gaia_dr3_ls_dr9_G_14_18.fits', overwrite=True)


cat = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/gaia_dr3_ls_dr9_G_14_18.fits'))
mask = cat['PHOT_BP_MEAN_MAG']-cat['PHOT_RP_MEAN_MAG']<2.4
cat = cat[mask]
SynthGaia = -0.509 + 0.459*cat['FLUX_G'] + 0.433*cat['FLUX_R'] + 0.163*cat['FLUX_Z']
G = 22.5-2.5*np.log10(SynthGaia)
cat['diff'] = cat['PHOT_G_MEAN_MAG'] - G


def get_stats_in_pixel(pix_idx):

    pix_list = pix_unique[pix_idx]

    hp_table = Table()
    hp_table['HPXPIXEL'] = pix_list
    # hp_table['RA'], hp_table['DEC'] = hp.pixelfunc.pix2ang(nside, pix_list, nest=False, lonlat=True)

    hp_table['mean'] = 0.
    hp_table['median'] = 0.
    hp_table['n_star'] = 0

    for index in np.arange(len(pix_idx)):
        idx = pixorder[pixcnts[pix_idx[index]]:pixcnts[pix_idx[index]+1]]
        v = cat['diff'][idx].copy()
        hp_table['mean'][index] = np.mean(v)
        hp_table['median'][index] = np.median(v)
        hp_table['n_star'][index] = len(v)

    return hp_table


nside = 128

pix_allobj = hp.pixelfunc.ang2pix(nside, cat['RA'], cat['DEC'], lonlat=True)
pix_unique, pixcnts = np.unique(pix_allobj, return_counts=True)
pixcnts = np.insert(pixcnts, 0, 0)
pixcnts = np.cumsum(pixcnts)
pixorder = np.argsort(pix_allobj)
n_process = 32
pix_idx_split = np.array_split(np.arange(len(pix_unique)), n_process)
with Pool(processes=n_process) as pool:
    res = pool.map(get_stats_in_pixel, pix_idx_split)

maps = vstack(res)
maps.sort('HPXPIXEL')

maps.write('/pscratch/sd/r/rongpu/tmp/gaia_dr3_ls_dr9_G_14_18_map_{}.fits'.format(nside), overwrite=True)
