from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

run_list_bad = [426, 431, 518, 553, 592, 703, 706, 736] 
run_list_replacement = [425, 430, 519, 552, 591, 704, 707, 737]

skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8.fits')

for run_bad, run_replacement in zip(run_list_bad, run_list_replacement):

    mask = skyrun['run']==run_bad
    band = skyrun['filter'][mask][0]

    bad_fn = os.path.join('sky_template_{}_{}.fits.fz'.format(band, run_bad))
    replacement_fn = os.path.join('sky_template_{}_{}.fits.fz'.format(band, run_replacement))

    print('cp {} {}'.format(replacement_fn, bad_fn))
