import sys
import numpy as np
from pathlib import Path

# sys.path.insert(0, "/home/jon/temp/codex_play")
# from stack_desi_restframe import make_wave_grid, stack_spectra

tids = [39089837617205054,39089837617196690,39089837617196796,39089837617196812]
tids = np.asarray(tids, dtype=np.int64)

# Fast TARGETID -> catalog row lookup
cat_targetids = np.asarray(cat["TARGETID"], dtype=np.int64)
row_by_tid = {int(tid): i for i, tid in enumerate(cat_targetids)}

spectra_paths = set()
z_by_targetid = {}

for tid in tids:
    tid = int(tid)

    if tid not in row_by_tid:
        print(f"Missing TARGETID in catalog: {tid}")
        continue

    index = row_by_tid[tid]
    z = float(cat["Z_cnn_rr"][index])

    if not np.isfinite(z) or z <= -1.0:
        print(f"Skipping invalid redshift for TARGETID {tid}: z={z}")
        continue

    uniqpix = int(cat["UNIQPIX"][index])
    pixgroup = int(uniqpix / 100)

    coadd_fn = (
        "/global/cfs/cdirs/desi/spectro/redux/matterhorn/spectra/special/other/"
        f"{pixgroup}/{uniqpix}/coadd-special-other-{uniqpix}.fits"
    )

    spectra_paths.add(coadd_fn)
    z_by_targetid[tid] = z

spectra_paths = sorted(spectra_paths)

rest_wave = make_wave_grid(800.0, 10000.0, 0.1)

stacked_flux, ncontrib, nspec_used = stack_spectra(
    spectra_paths=spectra_paths,
    redshifts=z_by_targetid,
    rest_wave_grid=rest_wave,
    cameras=None,
)

print(f"Stacked {nspec_used} spectra")

import matplotlib.pyplot as plt

plt.figure(figsize=(12, 4))
plt.plot(rest_wave, stacked_flux, lw=0.7)
plt.ylim(-1, 1)
plt.xlabel("Rest-frame wavelength [Angstrom]")
plt.ylabel("Mean flux")
plt.title(f"Rest-frame mean stack, N={nspec_used}")
plt.show()

