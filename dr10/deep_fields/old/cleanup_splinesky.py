from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool


n_processes = 4

# Load CCD list
ccd = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/dr10dev/deep_fields/survey-ccds-dr10-deep-fields-v1-defringed.fits'))
print(len(ccd))

_, idx = np.unique(ccd['expnum'], return_index=True)
exp = ccd[idx].copy()
print(len(exp))

expnum_list = [229707, 229714, 231597, 231606, 239013, 239014, 239023, 239039, 242475, 243877, 251149, 253162, 253359, 261959, 262227, 266864, 271191, 273346, 273637, 276026, 276028, 328717, 358528, 358533, 358536, 361622, 361648, 363793, 367215, 367234, 370185, 374572, 376728, 377467, 381579, 386030, 390479, 394219, 395520, 396626, 404467, 404775, 469194, 472394, 473993, 473994, 486857, 492478, 494676, 495331, 497358, 500271, 500309, 500495, 500508, 500776, 500780, 506656, 506796, 508822, 511338, 513310, 565129, 567450, 569167, 569618, 571484, 574285, 576152, 591428, 591447, 591493, 591763, 591770, 593164, 593746, 595070, 596471, 596799, 598252, 598981, 598989, 603993, 606301, 609593, 610006, 617097, 670052, 672370, 675283, 677333, 677334, 677340, 677341, 680560, 680577, 686880, 686886, 688493, 692781, 696013, 697623, 698411, 698418, 700332, 700380, 702073, 703224, 704142, 704150, 704786, 705145, 705185, 706142, 706149, 708216, 710924, 785421, 785439, 785813, 795315, 804768, 905991, 906038, 906074, 908580, 912507, 686923, 617649, 253196, 785337, 593167, 276397, 384194, 405171, 262252, 396348, 703218, 257574, 253421, 257581, 714833, 377558]
mask = np.in1d(exp['expnum'], expnum_list)
exp = exp[mask]
print(len(exp))


def check_splinesky(index):
    fn = '/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/calib/sky/' + exp['image_filename'][index].replace('.fits.fz', '-splinesky.fits')
    fn_dest = os.path.join('/pscratch/sd/r/rongpu/tractor/deep_fields/cosmos/attic/bad_splinesky/', os.path.basename(fn))
    cmd = 'mv {} {}'.format(fn, fn_dest)
    print(cmd)
    os.system(cmd)


print('Start!')
time_start = time.time()

with Pool(processes=n_processes) as pool:
    res = pool.map(check_splinesky, np.arange(len(exp)))

print('Done!', time.strftime('%H:%M:%S', time.gmtime(time.time() - time_start)))
