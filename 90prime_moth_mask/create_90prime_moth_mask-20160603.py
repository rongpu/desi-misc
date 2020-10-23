from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np


expnum_list = [75430048, 75430049, 75430050, 75430051, 75430052, 75430053, 75430054, 75430055, 75430056]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3600:, :1000] = True
    data['CCD3'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430057]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3350:, :1300] = True
    data['CCD3'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430073]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3950:, 1400:2350] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:1350, 250:2450] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430074]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3950:, 1400:2350] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:1350, 250:2450] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430075]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3950:, 1400:2350] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:1350, 350:2450] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430076]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3900:, 1400:2550] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:1300, 450:2500] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430077]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3900:, 1400:2550] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:1250, 450:2550] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430078]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3900:, 1400:2650] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:1200, 500:2650] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430079]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3800:, 1400:2750] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:1100, 700:2750] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430080]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3800:, 1400:2900] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:900, 750:2900] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430082]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3800:, 1400:3000] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:850, 900:3000] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430083]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3800:, 1400:3000] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:800, 900:3100] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430084]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[3800:, 1400:3000] = True
    data['CCD2'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[:800, 900:3000] = True
    data['CCD4'] = mask
    
    np.savez_compressed(output_path, **data)

expnum_list = [75430087]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[2500:3500, 3800:] = True
    data['CCD1'] = mask

    mask = np.zeros((4096, 4032), dtype=bool)
    mask[4096//2:3700, :4032//2] = True
    data['CCD2'] = mask

    np.savez_compressed(output_path, **data)

expnum_list = [75430088]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[4096//2:3700, :4032//2] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430089]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[4096//2:3400, :2250] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430090]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1900:3350, 250:2500] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430091]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1750:3100, 500:2800] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430092]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1700:3000, 700:3100] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430093]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1500:2800, 1100:3300] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430094]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1300:2700, 1350:3600] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430095]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1200:2600, 1400:3700] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430096]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[1250:2550, 1700:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430097]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[700:2450, 2000:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430098]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[350:2400, 2500:] = True
    mask[:3000, :350] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430099]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[350:2400, 2600:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430100]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[250:2250, 2800:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430101]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[250:2150, 2900:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430102]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[200:2150, 3200:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430103]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[200:1800, 3400:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430104]
for expnum in expnum_list:
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[200:1700, 3500:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430105, 75430106, 75430107, 75430108, 75430109]
for expnum in expnum_list:    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[300:1600, 3650:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)

expnum_list = [75430110, 75430111, 75430112, 75430113, 75430114, 75430115, 75430116, 75430117, 75430118, 75430119, 75430120, 75430121, 75430122, 75430123, 75430124, 75430125, 75430126]
for expnum in expnum_list:    
    output_path = '/global/u2/r/rongpu/temp/90prime_junk/mask_{}.npz'.format(expnum)
    data = {}
    mask = np.zeros((4096, 4032), dtype=bool)
    mask[350:1450, 3750:] = True
    data['CCD2'] = mask
    np.savez_compressed(output_path, **data)
