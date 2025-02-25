from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

##########################################################################################

f = open("/global/cfs/cdirs/desi/users/rongpu/plots/dr9dev/weird_psf/decam_largest_1st_moment.html", "w")
f.write('<html>\n')
f.write('<table>\n')

for ccd_index in ccd_index_list_plot:

    plot_fn1 = 'decam_largest_1st_moment/{}.png'.format(ccd['expnum'][ccd_index])

    f.write('<tr>\n')

    f.write('<td><a href=\'{}\'><img src=\'{}\' width=\'400\'></a></td>\n'.format(plot_fn1, plot_fn1))
    
    # f.write('<td><a href=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\'><img src=\'lg_regular_cp-sky_corrected/c4d_180221_072207_ooi_g_ls9.png\' width=\'400\'></a></td>\n')
    
    f.write('</tr>\n')
        
f.write('</table>\n')
f.close()
