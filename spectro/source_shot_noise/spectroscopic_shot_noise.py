import numpy as np

h = 6.62607015e-27   # erg s
c = 2.99792458e10    # cm/s

wavelength = 5000e-8  # cm (5000 Angstrom)
flux = 1e-17          # erg/s/cm^2/A, flat source spectrum
A_tel = 8.658e4       # cm^2; Mayall collecting area with DESI cage is 8.658 m2
                      # https://noirlab.edu/science/programs/kpno/telescopes/nicholas-mayall-4m-telescope/basic-parameters
tau = 0.30            # total system throughput at 5000 Angstrom
dlam = 0.8            # Angstrom/pixel, DESI extraction sampling
exptime = 1000.0      # s, nominal dark-time exposure
gain = 1.0            # e-/ADU

E_photon = h * c / wavelength
N_e = flux * dlam * A_tel * tau * exptime / E_photon
N_adu = N_e / gain
shot_noise_e = np.sqrt(N_e)
shot_noise_adu = shot_noise_e / gain

print(f"photon energy at 500 nm: {E_photon:.4e} erg")
print(f"detected electrons: {N_e:.2f} e-")
print(f"counts: {N_adu:.2f} ADU")
print(f"source shot noise: {shot_noise_adu:.2f} ADU")
