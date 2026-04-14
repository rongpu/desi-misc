import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp

dirname = '/pscratch/sd/r/rongpu/desi2/spectro/lae_redshift_efficiency_vs_flux'

lae_spectrum = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/raichoor/laelbg/lya-flim/desi-lya-flim-coadd-skyflux.fits', ext='LYAPROF'))
mwave, mflux = np.array(lae_spectrum['WAVE']), np.array(lae_spectrum['FLUX'])
plt.plot(mwave, mflux)
mask = (mwave>1210) & (mwave<1225)
mwave = mwave[mask]
mflux = mflux[mask]
plt.plot(mwave, mflux)
plt.show()
plt.show()

from scipy.integrate import trapezoid
print('Total flux: {:.1f} [10^-17 erg s^-1 cm^-2]'.format(trapezoid(mflux, mwave)))

fn_original = '/global/cfs/cdirs/desi/spectro/redux/loa/tiles/cumulative/1263/20240212/coadd-0-1263-thru20240212.fits'
b_wavelength = fitsio.read(fn_original, ext='B_WAVELENGTH')

for total_flux in [0.5, 1.0, 2.0]:

    plt.figure(figsize=(7, 6.5))
    for band_index, band in enumerate(['M411', 'M438', 'M464', 'M490', 'M517']):

        n_stack_arr = np.arange(1, 8)
        median_sigma_detect_arr = np.full(len(n_stack_arr), np.nan)

        for index, n_stack in enumerate(n_stack_arr):

            cat_stack = []
            for obs_index in [0, 1, 2, 3]:

                coadd_fn = os.path.join(dirname, 'coadd-lae-{:.1f}-{}-{}-{}-br_only.fits'.format(total_flux, obs_index, n_stack, band.lower()))
                redrock_fn = os.path.join(dirname, 'redrock-lae-{:.1f}-{}-{}-{}-br_only.fits'.format(total_flux, obs_index, n_stack, band.lower()))

                if not os.path.isfile(redrock_fn):
                    continue

                cat = Table(fitsio.read(redrock_fn, ext='REDSHIFTS'))
                tmp = Table(fitsio.read(redrock_fn, ext='FIBERMAP'))
                cat['Z_TRUE'] = tmp['TARGET_RA']
                cat['sigma_detect'] = -99.

                wave = fitsio.read(coadd_fn, ext='B_WAVELENGTH')
                flux_all = fitsio.read(coadd_fn, ext='B_FLUX')
                ivar_all = fitsio.read(coadd_fn, ext='B_IVAR')

                for ii in range(len(cat)):
                    redshift = cat['Z_TRUE'][index]
                    flux = flux_all[index]
                    ivar = ivar_all[index]
                    flux_norm_factor = 10.*total_flux / (redshift+1)
                    mflux1 = flux_norm_factor * mflux
                    mwave1 = mwave * (redshift+1)
                    b_flux_lae = np.interp(b_wavelength, mwave1, mflux1, left=0, right=0)
                    cat['sigma_detect'][ii] = np.sqrt(np.sum(ivar * b_flux_lae**2))

                cat_stack.append(cat)

            cat = vstack(cat_stack)

            median_sigma_detect_arr[index] = np.median(cat['sigma_detect'])

        plt.plot(n_stack_arr, median_sigma_detect_arr, '.-', label=band, color='C'+str(band_index))

    plt.legend(loc='lower right')
    plt.xlabel('N_exposure')
    plt.ylabel('Median detection sigmas')
    plt.grid(alpha=0.5)
    plt.title(r'Ly$\alpha$ flux: {:.1f}$\times 10^{{-16}} \ \mathrm{{erg}} \ \mathrm{{s}}^{{-1}} \ \mathrm{{cm}}^{{-2}}$'.format(total_flux))
    # plt.ylim(ymin=-0.01)
    plt.show()
