# Compute the initial fringe templates to get the initial fringe scales for normalization

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack, join
import fitsio

from pathlib import Path
from multiprocessing import Pool

from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter
from scipy import stats


time_start = time.time()

params = {'legend.fontsize': 'large',
          'axes.labelsize': 'large',
          'axes.titlesize': 'large',
          'xtick.labelsize': 'large',
          'ytick.labelsize': 'large',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)

fringe_dr9_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/fringe/DECam_CP-Fringe-Normed'
fringe_init_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/fringe_templates_init'
image_dir = '/global/cfs/cdirs/cosmo/staging/'
surveyccd_path = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/survey-ccds-z-fringe.fits'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/decam_ccd_blob_mask'

output_dir = '/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/frgscale_init'

# Load CCD list
ccd_columns = ['image_filename', 'image_hdu', 'ccdname', 'filter', 'expnum', 'exptime', 'ccdskycounts']
ccd = Table(fitsio.read(surveyccd_path, columns=ccd_columns))
print(len(ccd))

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx].copy()
str_loc = np.char.find(np.array(exp['image_filename'], dtype='str'), '/CP20')
exp['obs_date'] = np.array([exp['image_filename'][i][str_loc[i]+3:str_loc[i]+11] for i in range(len(exp))])
print('Total number of nights:', len(np.unique(exp['obs_date'])))
print('Number of images: {}'.format(len(exp)))


def check_frgscale(expnum):

    mask = exp['expnum']==expnum
    image_fn = exp['image_filename'][mask][0]
    frgscale_output_path = os.path.join(output_dir, image_fn.strip().replace('.fits.fz', '.txt'))

    if not os.path.isfile(frgscale_output_path):
        print('Missing!!!', frgscale_output_path)
        return None

    if os.stat(frgscale_output_path).st_size==0:
        print('Empty', frgscale_output_path)
        return None

    try:
        tmp = Table.read(frgscale_output_path, format='ascii.commented_header')
        print('expnum:', expnum, len(tmp))
    except:
        print('Error!!', frgscale_output_path)

    return tmp


print('=============================================================================')

expnum_list = np.array(exp['expnum'])

n_processes = 64
with Pool(processes=n_processes) as pool:
    res = pool.map(check_frgscale, expnum_list, chunksize=1)

# Remove None elements from the list
for index in range(len(res)-1, -1, -1):
    if res[index] is None:
        res.pop(index)

frgscales = vstack(res)
frgscales.write('/global/cfs/cdirs/desi/users/rongpu/dr10dev/fringe/data/frgscale_init/survey-ccds-dr10-deep-fields-v1-z-fringe-stacking-frgscales.fits')

print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))

