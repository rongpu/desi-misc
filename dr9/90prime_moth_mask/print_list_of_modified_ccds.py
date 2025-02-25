from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio


surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9m/survey-ccds-90prime-dr9.fits.gz'
ccd = Table(fitsio.read(surveyccd_path))

# Only keep unique exposures
_, idx = np.unique(ccd['expnum'], return_index=True)
ccd = ccd[idx]

########################################### 2016-01-11 ###########################################

expnum_min = 73990064
expnum_max = 73990184

mask = (ccd['expnum']>=expnum_min) & (ccd['expnum']<=expnum_max)
# print(np.sum(mask))
expnum_list = np.unique(ccd['expnum'][mask])
# print(len(expnum_list))
# print(expnum_list)

for expnum in expnum_list:
    print(expnum, 'CCD1')

########################################### 2016-06-03 ###########################################

expnum_list = [75430048, 75430049, 75430050, 75430051, 75430052, 75430053, 75430054, 75430055, 75430056]
for expnum in expnum_list:
    print(expnum, 'CCD3')

expnum_list = [75430057]
for expnum in expnum_list:
    print(expnum, 'CCD3')

expnum_list = [75430073]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430074]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430075]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430076]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430077]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430078]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430079]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430080]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430082]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430083]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430084]
for expnum in expnum_list:
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [75430087]
for expnum in expnum_list:
    print(expnum, 'CCD1')
    print(expnum, 'CCD2')


expnum_list = [75430088]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430089]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430090]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430091]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430092]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430093]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430094]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430095]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430096]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430097]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430098]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430099]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430100]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430101]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430102]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430103]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430104]
for expnum in expnum_list:
    print(expnum, 'CCD2')

expnum_list = [75430105, 75430106, 75430107, 75430108, 75430109]
for expnum in expnum_list:    
    print(expnum, 'CCD2')

expnum_list = [75430110, 75430111, 75430112, 75430113, 75430114, 75430115, 75430116, 75430117, 75430118, 75430119, 75430120, 75430121, 75430122, 75430123, 75430124, 75430125, 75430126]
for expnum in expnum_list:    
    print(expnum, 'CCD2')

########################################### 2017-06-19 ###########################################

expnum_list = [79240047, 79240048, 79240049, 79240050]
for expnum in expnum_list:
    print(expnum, 'CCD1')

expnum_list = [79240057, 79240058, 79240059, 79240060]
for expnum in expnum_list:
    print(expnum, 'CCD1')
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [79240061]
for expnum in expnum_list:
    print(expnum, 'CCD1')
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [79240062]
for expnum in expnum_list:
    print(expnum, 'CCD1')
    print(expnum, 'CCD2')
    print(expnum, 'CCD3')
    print(expnum, 'CCD4')

expnum_list = [79240063]
for expnum in expnum_list:
    print(expnum, 'CCD1')
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')

expnum_list = [79240064]
for expnum in expnum_list:
    print(expnum, 'CCD1')
    print(expnum, 'CCD2')
    print(expnum, 'CCD3')
    print(expnum, 'CCD4')

expnum_list = [79240065]
for expnum in expnum_list:
    print(expnum, 'CCD1')
    print(expnum, 'CCD2')
    print(expnum, 'CCD4')


