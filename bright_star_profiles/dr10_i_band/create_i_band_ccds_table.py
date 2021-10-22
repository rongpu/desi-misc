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

ccds = []
for fn in file_list:
    tt = Table(fitsio.read(fn))
    ccds.append(tt)

ccds = vstack(ccds)
print(len(ccds))
ccds.write('/global/cscratch1/sd/rongpu/temp/dr10_i_band_ccds.fits')
