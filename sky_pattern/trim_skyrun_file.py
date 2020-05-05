# Create a smaller list of exposure to speed up the initial blobmask run

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

############################ Exposures for creating templates ############################
##################################### Only 50 runs #######################################

n_run = 50
max_exposure = 50

skyrun = Table.read('/global/cscratch1/sd/schlafly/legacysurvey/skyrunsgoodcountexpnumv48dr8.fits')
print(len(skyrun))
mask = skyrun['ok']==True
skyrun = skyrun[mask]
print(len(skyrun))

idx_keep = []
for run in np.unique(skyrun['run']):
    idx = np.where(skyrun['run']==run)[0]
    if len(idx)>max_exposure:
        idx = idx[:max_exposure]
    idx_keep.append(idx)
idx_keep = np.concatenate(idx_keep)
skyrun = skyrun[idx_keep]
print(len(skyrun))

run_list = []
np.random.seed(321)
for band in ['g', 'r', 'z']:
    mask = skyrun['filter']==band
    run_list.append(np.random.choice(np.unique(skyrun['run'][mask]), size=n_run, replace=False))
run_list = np.concatenate(run_list)

mask = np.in1d(skyrun['run'], run_list)
skyrun = skyrun[mask]
print(len(skyrun))

for band in ['g', 'r', 'z']:
    mask = skyrun['filter']==band
    print(band, np.sum(mask), len(np.unique(skyrun['run'][mask])))
    a, b = np.unique(skyrun['run'][mask], return_counts=True)
    print(b.max(), b.min())

skyrun.write('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8_less.fits')

############################ Exposures for creating templates ############################
####################################### All runs #########################################

max_exposure = 50

skyrun = Table.read('/global/cscratch1/sd/schlafly/legacysurvey/skyrunsgoodcountexpnumv48dr8.fits')
print(len(skyrun))
mask = skyrun['ok']==True
skyrun = skyrun[mask]
print(len(skyrun))

idx_keep = []
for run in np.unique(skyrun['run']):
    idx = np.where(skyrun['run']==run)[0]
    if len(idx)>max_exposure:
        idx = idx[:max_exposure]
    idx_keep.append(idx)
idx_keep = np.concatenate(idx_keep)
skyrun = skyrun[idx_keep]
print(len(skyrun))

for band in ['g', 'r', 'z']:
    mask = skyrun['filter']==band
    print(band, np.sum(mask), len(np.unique(skyrun['run'][mask])))
    a, b = np.unique(skyrun['run'][mask], return_counts=True)
    print(b.max(), b.min())

skyrun.write('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8_less_than_50_exposures.fits')

############################ Exposures for testing sky subtraction ############################

n_exposure_in_each_band = 300

skyrun = Table.read('/global/cscratch1/sd/schlafly/legacysurvey/skyrunsgoodcountexpnumv48dr8.fits')
print(len(skyrun))

idx_keep = []

np.random.seed(12)
for band in ['g', 'r', 'z']:
    idx = np.where(skyrun['filter']==band)[0]
    idx = np.random.choice(idx, size=n_exposure_in_each_band, replace=False)
    idx_keep.append(idx)
idx_keep = np.concatenate(idx_keep)
skyrun = skyrun[idx_keep]
print(len(skyrun))

for band in ['g', 'r', 'z']:
    mask = skyrun['filter']==band
    print(band, np.sum(mask), len(np.unique(skyrun['run'][mask])))

skyrun.write('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8_random_subset.fits')

