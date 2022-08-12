from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

n_subsets = 11

sub_dir_x = '/pscratch/sd/r/rongpu/tractor/deep_field_subsets/cosmos/sub-x'

for subset in range(n_subsets):

    sub_dir = '/pscratch/sd/r/rongpu/tractor/deep_field_subsets/cosmos/sub-{}'.format(subset)

    if os.path.isdir(sub_dir):
        continue

    os.makedirs(sub_dir)
    os.system('rsync -av {}/* {}/'.format(sub_dir_x, sub_dir))
    os.system('ln -s /global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/cosmos/survey-ccds-dr10-deep-fields-v1-defringed-subset-{}.fits {}/survey-ccds-dr10-deep-fields-v1-defringed-subset-{}.fits'.format(subset, sub_dir, subset))
    os.system('ln -s /global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_field_subsets/cosmos/survey-ccds-dr10-deep-fields-v1-defringed-subset-{}.kd.fits {}/survey-ccds-dr10-deep-fields-v1-defringed-subset-{}.kd.fits'.format(subset, sub_dir, subset))
