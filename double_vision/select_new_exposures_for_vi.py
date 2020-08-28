from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
# import matplotlib.pyplot as plt
import numpy as np

# DECam
expnum_list = glob.glob('/global/cscratch1/sd/rongpu/temp/double_vision/decam_*')
expnum_list = [int(os.path.basename(tmp).replace('decam_', '')) for tmp in expnum_list]
for expnum in expnum_list:
    os.system('cp /global/cscratch1/sd/rongpu/temp/double_vision_plots_round_2/image_decam_{}_* /global/cscratch1/sd/rongpu/temp/double_vision_plots_round_2_new'.format(expnum))

# 90prime
expnum_list = glob.glob('/global/cscratch1/sd/rongpu/temp/double_vision/90prime_*')
expnum_list = [int(os.path.basename(tmp).replace('90prime_', '')) for tmp in expnum_list]
for expnum in expnum_list:
    os.system('cp /global/cscratch1/sd/rongpu/temp/double_vision_plots_90prime_round_2/image_90prime_{}_* /global/cscratch1/sd/rongpu/temp/double_vision_plots_90prime_round_2_new'.format(expnum))

# mosaic
expnum_list = glob.glob('/global/cscratch1/sd/rongpu/temp/double_vision/mosaic_*')
expnum_list = [int(os.path.basename(tmp).replace('mosaic_', '')) for tmp in expnum_list]
for expnum in expnum_list:
    os.system('cp /global/cscratch1/sd/rongpu/temp/double_vision_plots_mosaic_round_2/image_mosaic_{}_* /global/cscratch1/sd/rongpu/temp/double_vision_plots_mosaic_round_2_new'.format(expnum))
