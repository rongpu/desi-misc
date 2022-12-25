# Perlmutter scratch only 8-weeks!

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

fns = glob.glob('/pscratch/sd/r/rongpu/dr10dev/fringe_corrected_images/decam/CP/*/*/*.fits.fz')
print(len(fns))

for fn in fns:
    tmp = fitsio.read_header(fn, ext=0)
