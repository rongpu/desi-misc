# Convert the files from .csv.gz to .fits

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
from multiprocessing import Pool

output_dir = '/pscratch/sd/r/rongpu/gaia_dr3/xp_continuous_mean_spectrum'
fns = sorted(glob.glob('/pscratch/sd/r/rongpu/gaia_dr3/xp_continuous_mean_spectrum/hardtouseformat/XpContinuousMeanSpectrum_*.csv.gz'))
print(len(fns))


def do_something(fn):

    output_fn = os.path.join(output_dir, os.path.basename(fn).replace('.csv.gz', '.fits'))
    if os.path.isfile(output_fn):
        print(output_fn, 'already exists!')
        return None

    try:
        cat = Table.read(fn, format='ascii.ecsv')
    except AssertionError:
        return None

    # convert int8 to int16 because int8 is not supported by FITS
    for col in cat.colnames:
        if cat[col].dtype==np.int8:
            cat[col] = np.array(cat[col]).astype(np.int16)

    # convert "object" type to numpy array so it can be saved in FITS
    for col in ['bp_coefficients', 'bp_coefficient_errors', 'bp_coefficient_correlations', 'rp_coefficients', 'rp_coefficient_errors', 'rp_coefficient_correlations']:
        cat[col] = np.array(list(cat[col]))

    cat.write(output_fn)

    gc.collect()

    return None


n_process = 16
with Pool(processes=n_process) as pool:
    res = pool.map(do_something, fns, chunksize=1)

