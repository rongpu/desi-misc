from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from multiprocessing import Pool

n_processes = 16

surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'

ccd_columns = ['image_filename', 'image_hdu', 'expnum', 'filter']
ccd = Table(fitsio.read(surveyccd_path, ccd_columns=ccd_columns))

# unique exposures only
ccd.sort('expnum')
mask = np.concatenate([[True], np.diff(ccd['expnum'])!=0])
ccd = ccd[mask]

print(len(ccd), 'exposures')

# for index in range(len(ccd)):
def check_integrity(index):

    if index%1000==0:
        print(index, '/', len(ccd))

    expnum = ccd['expnum'][index]
    blob_path = os.path.join(blob_dir, 'blob_mask', ccd['image_filename'][index].replace('.fits.fz', '-blobmask.npz').strip())
    
    if os.stat(blob_path).st_size == 0:
        # continue
        return None

    try:
        tmp = np.load(blob_path)
    except:
        print(expnum, blob_path)

def main():

    with Pool(processes=n_processes) as pool:
        res = pool.map(check_integrity, np.arange(len(ccd)))

    print('Done!!!!!!!!!!!!!!!!!!!!!')

if __name__=="__main__":
    main()

