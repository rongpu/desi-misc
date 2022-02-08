# Write images without the CP fringe correction

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
from astropy.io import fits
import fitsio

from multiprocessing import Pool


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
                  'N29': 60, 'N30': 61, 'N31': 62}
ccdnamenumdict_inv = {aa: bb for bb, aa in ccdnamenumdict.items()}

ccdnum_list = [ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
               18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
               35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
               52, 53, 54, 55, 56, 57, 58, 59, 60, 62]


image_dir = '/global/project/projectdirs/cosmo/staging'
output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/undo_cp_fringe_corr'
cp_fringe_dir = '/global/homes/d/djschleg/cosmo/staging/decam/DECam_CP-Fringe'

fringe_list_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/misc/fringe_list.txt'
with open(fringe_list_path, "r") as f:
    fringe_list = f.read().splitlines()

# Load old fringe image
fringe_old_dict = {}
for ccdnum in ccdnum_list:
    # skip N30 and S7
    if ccdnum==61 or ccdnum==31:
        continue
    fringe_old_path = os.path.join(cp_fringe_dir, 'DES17B_20180103_908c062-z-{}_frg.fits'.format(str(ccdnum).zfill(2)))
    fringe_old = fits.getdata(fringe_old_path)
    # remove the edge pixels
    fringe_old = fringe_old[1:4095, 1:2047]
    fringe_old_dict[ccdnum] = fringe_old.copy()


def write_new_image(fn):

    img_fn = os.path.join(image_dir, fn)
    img_fn_write = os.path.join(output_dir, fn)

    print(img_fn_write)

    if not os.path.exists(os.path.dirname(img_fn_write)):
        try:
            os.makedirs(os.path.dirname(img_fn_write))
        except:
            pass

    hdul_r = fitsio.FITS(img_fn, mode='r')
    hdul_w = fitsio.FITS(img_fn_write, mode='rw', clobber=True)

    for hdu_index in range(len(hdul_r)):
        if hdu_index==0:
            hdr = hdul_r[hdu_index].read_header()
            hdul_w.write(data=None, header=hdr)
        else:
            hdr = hdul_r[hdu_index].read_header()
            img = hdul_r[hdu_index].read()
            ccdname = hdul_r[hdu_index].get_extname()
            ccdnum = ccdnamenumdict[ccdname.strip()]

            # Back out the exisiting fringe correction
            if 'FRGSCALE' in hdr:  # for some CCD(s) (S7?) no FRGSCALE exist in the original header, and no correct is done here
                frgscale_old = hdr['FRGSCALE']
                # print(ccdname, frgscale_old)
                fringe_old = fringe_old_dict[ccdnum]
                img += fringe_old * frgscale_old
                hdr.delete('FRGSCALE')
                hdr.delete('FRINGE')
            else:
                print(ccdname, 'no fringe')

            hdr.delete('CHECKSUM')
            hdr.delete('DATASUM')

            hdul_w.write(data=img, header=hdr, extname=ccdname, compress='rice')

    hdul_r.close()
    hdul_w.close()


n_processes = 64
with Pool(processes=n_processes) as pool:
    res = pool.map(write_new_image, fringe_list)


