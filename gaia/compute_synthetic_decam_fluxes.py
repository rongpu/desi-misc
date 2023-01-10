# Convert the files from .csv.gz to .fits

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool
from gaiaxpy import generate, PhotometricSystem

output_dir = '/global/cfs/cdirs/desi/users/rongpu/data/gaia_dr3/xp_synthetic_decam_photometry'
fns = sorted(glob.glob('/pscratch/sd/r/rongpu/gaia_dr3/xp_continuous_mean_spectrum/XpContinuousMeanSpectrum_*.fits'))
print(len(fns))


def do_something(fn):

    output_fn = os.path.join(output_dir, os.path.basename(fn).replace('XpContinuousMeanSpectrum_', 'XpSyntheticDECam_'))
    if os.path.isfile(output_fn):
        print(output_fn, 'already exists!')
        return None

    cat = Table(fitsio.read(fn))

    # Workaround to make the multidimensional arrays pandas-compatible
    cat = vstack([cat, cat[:1].copy()])
    for col in ['bp_coefficients', 'bp_coefficient_errors', 'bp_coefficient_correlations', 'rp_coefficients', 'rp_coefficient_errors', 'rp_coefficient_correlations']:
        tmp = list(np.array(cat[col]))
        tmp.pop()
        tmp += [np.array([0])]
        cat[col] = tmp
    cat = cat[:-1]
    print(len(cat))

    cat = cat.to_pandas()

    phot_system = PhotometricSystem.DECam
    photom = generate(cat, photometric_system=phot_system, save_file=False)
    photom = Table.from_pandas(photom)

    print(np.allclose(photom['Decam_mag_g'], (-2.5*np.log10(photom['Decam_flux_g']))-56.1),
          np.allclose(photom['Decam_mag_r'], (-2.5*np.log10(photom['Decam_flux_r']))-56.1),
          np.allclose(photom['Decam_mag_i'], (-2.5*np.log10(photom['Decam_flux_i']))-56.1),
          np.allclose(photom['Decam_mag_z'], (-2.5*np.log10(photom['Decam_flux_z']))-56.1),
          np.allclose(photom['Decam_mag_Y'], (-2.5*np.log10(photom['Decam_flux_Y']))-56.1))

    for col in ['Decam_flux_g', 'Decam_flux_r', 'Decam_flux_i', 'Decam_flux_z', 'Decam_flux_Y', 'Decam_flux_error_g', 'Decam_flux_error_r', 'Decam_flux_error_i', 'Decam_flux_error_z', 'Decam_flux_error_Y']:
        photom[col.replace('Decam_', '')] = photom[col] * 10**31.44

    print(np.allclose(photom['Decam_mag_g'], (22.5-2.5*np.log10(photom['flux_g']))),
          np.allclose(photom['Decam_mag_r'], (22.5-2.5*np.log10(photom['flux_r']))),
          np.allclose(photom['Decam_mag_i'], (22.5-2.5*np.log10(photom['flux_i']))),
          np.allclose(photom['Decam_mag_z'], (22.5-2.5*np.log10(photom['flux_z']))),
          np.allclose(photom['Decam_mag_Y'], (22.5-2.5*np.log10(photom['flux_Y']))))

    photom = photom[['source_id', 'flux_g', 'flux_r', 'flux_i', 'flux_z', 'flux_Y', 'flux_error_g', 'flux_error_r', 'flux_error_i', 'flux_error_z', 'flux_error_Y']]
    photom.write(output_fn)

    return None


n_process = 16
with Pool(processes=n_process) as pool:
    res = pool.map(do_something, fns, chunksize=1)


