from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits
from astropy import wcs

input_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/data'
output_dir = '/global/project/projectdirs/desi/users/rongpu/dr9/fringe/DECam_CP-Fringe'

# http://www.ctio.noao.edu/noao/sites/default/files/DECam/DECamOrientation.png
ccdnamenumdict = {'S1': 25, 'S2': 26, 'S3': 27, 'S4':28,
                  'S5': 29, 'S6': 30, 'S7': 31,
                  'S8': 19, 'S9': 20, 'S10': 21, 'S11': 22, 'S12': 23,
                  'S13': 24,
                  'S14': 13, 'S15': 14, 'S16': 15, 'S17': 16, 'S18': 17,
                  'S19': 18,
                  'S20': 8, 'S21': 9, 'S22': 10, 'S23': 11, 'S24': 12,
                  'S25': 4, 'S26': 5, 'S27': 6, 'S28': 7,
                  'S29': 1, 'S30': 2, 'S31': 3,
                  'N1': 32, 'N2': 33, 'N3': 34, 'N4': 35,
                  'N5': 36, 'N6': 37, 'N7': 38,
                  'N8': 39, 'N9': 40, 'N10': 41, 'N11': 42, 'N12': 43,
                  'N13': 44,
                  'N14': 45, 'N15': 46, 'N16': 47, 'N17': 48, 'N18': 49,
                  'N19': 50,
                  'N20': 51, 'N21': 52, 'N22': 53, 'N23': 54, 'N24': 55,
                  'N25': 56, 'N26': 57, 'N27': 58, 'N28': 59,
                  'N29': 60, 'N30': 61, 'N31': 62,
                  }

# HDU index (image_hdu in survey-ccd table) to CCD name
hdu2ccdname = {1: 'S29', 2: 'S30', 3: 'S31', 4: 'S25', 5: 'S26', 6: 'S27', 7: 'S28', 
                8: 'S20', 9: 'S21', 10: 'S22', 11: 'S23', 12: 'S24', 13: 'S14', 14: 'S15', 
                15: 'S16', 16: 'S17', 17: 'S18', 18: 'S19', 19: 'S8', 20: 'S9', 21: 'S10', 
                22: 'S11', 23: 'S12', 24: 'S13', 25: 'S1', 26: 'S2', 27: 'S3', 28: 'S4', 
                29: 'S5', 30: 'S6', 31: 'S7', 32: 'N1', 33: 'N2', 34: 'N3', 35: 'N4', 36: 'N5', 
                37: 'N6', 38: 'N7', 39: 'N8', 40: 'N9', 41: 'N10', 42: 'N11', 43: 'N12', 
                44: 'N13', 45: 'N14', 46: 'N15', 47: 'N16', 48: 'N17', 49: 'N18', 50: 'N19', 
                51: 'N20', 52: 'N21', 53: 'N22', 54: 'N23', 55: 'N24', 56: 'N25', 57: 'N26', 
                58: 'N27', 59: 'N28', 60: 'N29', 61: 'N31'}

for ii, hdu_index in enumerate(range(1, 62)):

    if hdu_index==31:
        print('skipping S7')
        # sys.exit()
        continue

    print(hdu_index)

    data = np.load(os.path.join(input_dir, 'fringe_smooth_{}.npy'.format(hdu_index)))

    data_padded = np.zeros((data.shape[0]+2, data.shape[1]+2))
    data_padded[1:data_padded.shape[0]-1, 1:data_padded.shape[1]-1] = data
    
    mask = ~np.isfinite(data_padded)
    data_padded[mask] = 0
    data_padded = data_padded.astype(np.float32)

    hdu = fits.PrimaryHDU(data_padded)
    hdul = fits.HDUList([hdu])
    hdul.writeto(os.path.join(output_dir, 'DECam_z_frg_{}_{}_hdu{}.fits'.format(hdu2ccdname[hdu_index], str(ccdnamenumdict[hdu2ccdname[hdu_index]]).zfill(2), str(hdu_index).zfill(2))))

