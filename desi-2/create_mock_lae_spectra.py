# Based on /Users/rongpu/git/desi-misc/spectro/sky_fiber_noise/create_sky_fiber_coadds_with_lae.py

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

ibis_wave = {'M411': [4007, 4261.5], 'M438': [4261.5, 4522.0], 'M464': [4522.0, 4777.0], 'M490': [4777.0, 5039.5], 'M517': [5039.5, 5298]}

# spectro_config = 'b_only'
spectro_config = 'br_only'

output_dir = '/pscratch/sd/r/rongpu/desi2/spectro/lae_redshift_efficiency_vs_flux'

spec_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/spectro/sky_spectra/sky_spectra_different_nights_202401_202403.fits'

lae_spectrum = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/raichoor/laelbg/lya-flim/desi-lya-flim-coadd-skyflux.fits', ext='LYAPROF'))
wave, flux = np.array(lae_spectrum['WAVE']), np.array(lae_spectrum['FLUX'])

from scipy.integrate import trapezoid
print('Total flux: {:.1f} [10^-17 erg s^-1 cm^-2]'.format(trapezoid(flux, wave)))

# ibis_center_wave ={'M411': 4144, 'M438': 4396, 'M464': 4653, 'M490': 4911, 'M517': 5169}
# lae_ibis_mag = {'M411': 26.718704415635052, 'M438': 26.567682531084806, 'M464': 26.35383567788979, 'M490': 26.21582126791786, 'M517': 26.019661176572964}
# mag_norm = 25.0

fn_original = '/global/cfs/cdirs/desi/spectro/redux/loa/tiles/cumulative/1263/20240212/coadd-0-1263-thru20240212.fits'
b_wavelength = fitsio.read(fn_original, ext='B_WAVELENGTH')


specall = Table(fitsio.read(spec_fn))
print(len(specall))
fibers = np.unique(specall['FIBER'])
print(len(fibers))

# random shuffle
np.random.seed(777)
idx = np.random.choice(len(specall), size=len(specall), replace=False)
specall = specall[idx]
print(len(specall))

# Trim to 28 exposures per fiber and assign obs_index
n_obs = 4  # number of independent observations per fiber
n_stack_max = 7
n_exp_max = n_obs * n_stack_max
specall_stack = []
for fiber in fibers:
    mask = specall['FIBER']==fiber
    tmp = specall[mask][:n_exp_max].copy()
    tmp['obs_index'] = np.repeat(np.arange(n_obs), n_stack_max)
    specall_stack.append(tmp)
specall = vstack(specall_stack, join_type='exact')
print(len(specall))

# Use the same set of redshifts for all n_stack's
# np.random.seed(888)  # for obs_index=0
# np.random.seed(123)  # for obs_index=1
# np.random.seed(666)  # for obs_index=2
np.random.seed(100)  # for obs_index=3
redshifts_dict = {}
for band in ['M411', 'M438', 'M464', 'M490', 'M517']:
    lae_wave = np.random.uniform(ibis_wave[band][0], ibis_wave[band][1], 500)
    redshifts_dict[band] = lae_wave/1215.67-1

for total_flux in [1.0, 0.5, 2.0]:  # 10^-16 erg/s/cm2
# for total_flux in [1.0]:  # 10^-16 erg/s/cm2

    for band in ['M411', 'M438', 'M464', 'M490', 'M517']:

        redshifts = redshifts_dict[band]

        b_flux_lae_all = []
        for redshift in redshifts:
            wave1 = wave * (redshift+1)
            flux_norm_factor = 10.*total_flux / (redshift+1)
            flux1 = flux_norm_factor * flux
            b_flux_lae = np.interp(b_wavelength, wave1, flux1, left=0, right=0)
            # print('z = {:.5f}  Total flux: {:.1f} [10^-17 erg s^-1 cm^-2]'.format(redshift, trapezoid(b_flux_lae, b_wavelength)))
            b_flux_lae_all.append(b_flux_lae)
        b_flux_lae_all = np.array(b_flux_lae_all)

        # for obs_index in range(4):
        for obs_index in [3]:
            
            for n_stack in np.arange(1, 8):
            # for n_stack in [1]:
                fn_output = os.path.join(output_dir, 'coadd-lae-{:.1f}-{}-{}-{}.fits'.format(total_flux, obs_index, n_stack, band.lower()))

                if spectro_config=='b_only':
                    fn_output = fn_output.replace('.fits', '-b_only.fits')
                elif spectro_config=='br_only':
                    fn_output = fn_output.replace('.fits', '-br_only.fits')

                if os.path.isfile(fn_output):
                    print('file already exists; skip')
                    print(fn_output)
                    continue

                spec_stack = []
                for fiber in fibers:
                    spec = Table()
                    spec['FIBER'] = [fiber]
                    mask = (specall['obs_index']==obs_index) & (specall['FIBER']==fiber)
                    spec1 = specall[mask]
                    for camera in ['B', 'R', 'Z']:
                        ffivar = np.sum(spec1[camera+'_IVAR'][:n_stack] * spec1[camera+'_MSK'][:n_stack], axis=0)
                        ff = np.sum(spec1[camera+'_FLUX'][:n_stack] * spec1[camera+'_IVAR'][:n_stack] * spec1[camera+'_MSK'][:n_stack], axis=0) / ffivar
                        ff[np.isnan(ff)] = 0.
                        msk = np.sum(np.array(spec1[camera+'_MSK'][:n_stack]), axis=0)>0
                        spec[camera+'_FLUX'] = [ff]
                        spec[camera+'_IVAR'] = [ffivar]
                        spec[camera+'_MSK'] = [msk]
                        spec[camera+'_RESOLUTION'] = [spec1[camera+'_RESOLUTION'][0]]
                    spec_stack.append(spec)
                spec = vstack(spec_stack, join_type='exact')

                fits = fitsio.FITS(fn_original)

                with fitsio.FITS(fn_output, 'rw') as f:
                    for index in range(len(fits)):
                        extname = fits[index].get_extname()
                        data, header = fitsio.read(fn_original, ext=index, header=True)
                        if extname in ['FIBERMAP', 'EXP_FIBERMAP', 'SCORES']:
                            for col in data.dtype.names:
                                data[col] = data[col][0]
                            data['TARGETID'] = np.arange(500)
                            if extname=='FIBERMAP':
                                data['FIBER'] = np.array(spec['FIBER'])
                                data['TARGET_RA'] = redshifts  # use the TARGET_RA column to store redshifts
                        elif '_MASK' in extname:
                            data = np.array(~spec[extname[0]+'_MSK']).astype(np.int32) * 2**1  # ALLBADPIX
                        elif extname[1:] in ['_FLUX', '_IVAR', '_MSK', '_RESOLUTION']:
                            data = np.array(spec[extname])
                        if spectro_config=='b_only' and (extname.startswith('R_') or extname.startswith('Z_')):  # skip R and Z cameras
                            continue
                        if spectro_config=='br_only' and extname.startswith('Z_'):  # skip Z camera
                            continue
                        if extname=='B_WAVELENGTH':
                            assert np.all(data==b_wavelength)
                        if extname=='B_FLUX':
                            data += b_flux_lae_all
                        f.write(data, extname=extname, header=header)

                fits.close()
