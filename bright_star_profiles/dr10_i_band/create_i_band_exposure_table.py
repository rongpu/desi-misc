from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

file_list = glob.glob('/global/cscratch1/sd/dstn/dr10pre/zpt/decam/DECam_CP-DR10c/*/*_ooi_i*annotated.fits')
print(len(file_list))

exp = []
for fn in file_list:
    tt = Table(fitsio.read(fn))
    mask = tt['ccdname']=='S10'
    mask &= tt['filter']=='i'
    if np.sum(mask)==1:
        tt = tt[mask]
        exp.append(tt)
    elif np.sum(mask)>0:
        raise ValueError

exp = vstack(exp)
exp.write('/global/cscratch1/sd/rongpu/temp/dr10_i_band_exposures.fits')
