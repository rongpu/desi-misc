from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np

expnum_list = [79240047, 79240048, 79240049, 79240050]

for expnum in expnum_list:
    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)

    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:4096//2, :4032//2] = True
    data['CCD1'] = mask

    np.savez_compressed(output_path, **data)

    
expnum_list = [79240057, 79240058, 79240059, 79240060]

for expnum in expnum_list:    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)

    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1800:3000, 2016:] = True
    data['CCD1'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3000:, 500:3000] = True
    mask[1800:2900, 3600:] = True
    mask[3900:, 3600:] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:250, 200:2100] = True
    mask[:1100, 3500:] = True
    data['CCD4'] = mask

    np.savez_compressed(output_path, **data)


expnum_list = [79240061]

for expnum in expnum_list:    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)

    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1800:3000, 2016:] = True
    data['CCD1'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3000:, 500:3000] = True
    mask[2500:, 2500:] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:250, 200:2100] = True
    mask[:1100, 3500:] = True
    data['CCD4'] = mask

    np.savez_compressed(output_path, **data)


expnum_list = [79240062]

for expnum in expnum_list:    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)

    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1800:3000, 2016:] = True
    data['CCD1'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3000:, 500:3000] = True
    mask[3900:, 3600:] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1800:3600, :1000] = True
    data['CCD3'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:250, 200:2100] = True
    mask[:1100, 3500:] = True
    data['CCD4'] = mask

    np.savez_compressed(output_path, **data)


expnum_list = [79240063]

for expnum in expnum_list:    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)

    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1800:3000, 2016:] = True
    data['CCD1'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3000:, 200:3500] = True
    mask[3800:, 3500:] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:250, 200:2100] = True
    mask[:3000, 3000:] = True
    data['CCD4'] = mask

    np.savez_compressed(output_path, **data)


expnum_list = [79240064]

for expnum in expnum_list:    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)

    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1800:3000, 2016:] = True
    data['CCD1'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[2900:, 500:3000] = True
    mask[3900:, 3600:] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1000:3000, 2000:] = True
    data['CCD3'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:250, 200:2100] = True
    mask[:1100, 3500:] = True
    data['CCD4'] = mask

    np.savez_compressed(output_path, **data)


expnum_list = [79240065]

for expnum in expnum_list:    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)

    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1000:3000, 2000:] = True
    data['CCD1'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[2900:, 500:3000] = True
    mask[3900:, 3600:] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:250, 200:2100] = True
    mask[:1100, 3500:] = True
    data['CCD4'] = mask

    np.savez_compressed(output_path, **data)