import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord

hsc = Table(fitsio.read('/global/cfs/cdirs/desi/target/analysis/truth/parent/hsc-pdr3-dud-rev-no_mag_limit-reduced.fits',
                        columns=['ra', 'dec', 'g_cmodel_flux', 'r_cmodel_flux', 'i_cmodel_flux', 'z_cmodel_flux', 'y_cmodel_flux', 'a_g', 'a_r', 'a_i', 'a_z', 'a_y']))
print(len(hsc))

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    hsc['gmag_hsc_dereddened'] = -2.5 * np.log10(hsc['g_cmodel_flux']/3630.78) + 22.5 - hsc['a_g']
    hsc['rmag_hsc_dereddened'] = -2.5 * np.log10(hsc['r_cmodel_flux']/3630.78) + 22.5 - hsc['a_r']
    hsc['imag_hsc_dereddened'] = -2.5 * np.log10(hsc['i_cmodel_flux']/3630.78) + 22.5 - hsc['a_i']
    hsc['zmag_hsc_dereddened'] = -2.5 * np.log10(hsc['z_cmodel_flux']/3630.78) + 22.5 - hsc['a_z']
    hsc['ymag_hsc_dereddened'] = -2.5 * np.log10(hsc['y_cmodel_flux']/3630.78) + 22.5 - hsc['a_y']


for fn in [
'/global/cfs/cdirs/desi/users/raichoor/laelbg/clauds/dr/COSMOS_11bands-SExtractor-Lephare.fits',
'/global/cfs/cdirs/desi/users/raichoor/laelbg/clauds/dr/XMMLSS_11bands-SExtractor-Lephare.fits',
'/global/cfs/cdirs/desi/users/raichoor/laelbg/clauds/dr-offset/COSMOS_11bands-SExtractor-Lephare-offset.fits',
'/global/cfs/cdirs/desi/users/raichoor/laelbg/clauds/dr-offset/XMMLSS_11bands-SExtractor-Lephare-offset.fits',
]:

    cat = Table.read(fn)

    ramin, ramax, decmin, decmax = cat['RA'].min()-0.01, cat['RA'].max()+0.01, cat['DEC'].min()-0.01, cat['DEC'].max()+0.01
    mask = (hsc['ra']>ramin) & (hsc['ra']<ramax) & (hsc['dec']>decmin) & (hsc['dec']<decmax)
    hsc1 = hsc[mask].copy()
    print(len(hsc1))

    idx1, idx2, d2d, d_ra, d_dec = match_coord(hsc1['ra'], hsc1['dec'], cat['RA'], cat['DEC'], search_radius=1., plot_q=True)
    print(len(idx1)/len(cat))
    print(len(idx1)/len(hsc1))
    # hsc1 = hsc1[idx1]
    # cat = cat[idx2]

    for col in ['gmag_hsc_dereddened', 'rmag_hsc_dereddened', 'imag_hsc_dereddened', 'zmag_hsc_dereddened', 'ymag_hsc_dereddened']:
        cat[col] = np.nan
        cat[col][idx2] = hsc1[col][idx1]

    cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/misc/clauds/'+os.path.basename(fn).replace('.fits', '-HSC_PDR3_dud.fits'))


#########################################################################################################################################################

cat = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/misc/clauds/COSMOS_11bands-SExtractor-Lephare-HSC_PDR3_dud.fits'))

mask = np.isfinite(cat['rmag_hsc_dereddened'])
print(np.sum(mask), np.sum(mask)/len(mask))
mask &= cat['rmag_hsc_dereddened'] < 24.5
print(np.sum(mask), np.sum(mask)/len(mask))
mask &= cat['rmag_hsc_dereddened'] - cat['r'] < 0.5
print(np.sum(mask), np.sum(mask)/len(mask))

cat = cat[mask]
cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/misc/clauds/COSMOS_11bands-SExtractor-Lephare-HSC_PDR3_dud-clean.fits')
