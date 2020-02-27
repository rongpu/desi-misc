# Get the a list of science exposures for a specficic date and tile
# Print desi_coadd_spectra commands for single exposure coadds

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

# redux_dir = '/global/cfs/cdirs/desi/spectro/redux/daily'
# output_dir = '/global/cscratch1/sd/rongpu/desi/minisv2/single_exp_coadd'
redux_dir = os.getenv('REDUXDIR')
output_dir = os.getenv('OUTDIR')

obsdate_list = ['20200219']
overwrite = False

# tileid_list = None  # no restriction on tiles
tileid_list = [70003]
# tileid_list = [70003, 70004, 70005]

################################## Get list of exposures ##################################

exposure_dir_list = []
for obsdate in obsdate_list:
    exposure_dir_list += glob.glob(os.path.join(redux_dir, 'exposures', obsdate, '*'))

cframe_list = []

# Get a list of all science exposures
for exposure_dir in exposure_dir_list:

    cframes = glob.glob(os.path.join(exposure_dir, 'cframe-*'))

    if len(cframes)>0:

        if tileid_list is None:
            cframe_list += cframes
        else:
            # only need to check one cframe file in the exposure
            with fitsio.FITS(cframes[0]) as f:
                if f[0].read_header()['TILEID'] in tileid_list:
                    cframe_list += cframes

cframe_list = sorted(cframe_list)

################################## Print desi_coadd_spectra commands ##################################

print('\n################# desi_coadd_spectra commands: #################\n')

# list of exposures after the TILEID cut
exposure_dir_list = list(np.unique([os.path.dirname(cframe) for cframe in cframe_list]))

output_argument_list = []

for exposure_dir in exposure_dir_list:

    exposure = os.path.basename(exposure_dir)

    for petal_loc in range(10):
        
        cframes = glob.glob(os.path.join(exposure_dir, 'cframe-?{}-*').format(petal_loc))

        if len(cframes)==3:

            with fitsio.FITS(cframes[0]) as f:
                tileid = f[0].read_header()['TILEID']

            input_argument = os.path.join(exposure_dir.replace(redux_dir, '$REDUXDIR'), 'cframe-[brz]{}-*.fits').format(petal_loc)
            output_argument = os.path.join('$OUTDIR', str(tileid), 'coadd-{}-{}.fits'.format(petal_loc, exposure))
            output_argument_list.append(output_argument)

            if os.path.isfile(output_argument) and (not overwrite):
                print('\nWarninig: {} already exists!\n'.format(output_argument))
                continue

            print('time desi_coadd_spectra --coadd-cameras -i {} -o {}'.format(input_argument, output_argument))

        elif len(cframes)>0 and len(cframes)<3:

            print('\nWarning: less than three cframes files (ignored):\n', cframes, '\n')

################################## Print redrock commands ##################################

print('\n################# redrock commands: #################\n')

for output_argument in output_argument_list:

    rrdesi_argument_redrock = output_argument.replace('coadd', 'redrock').replace('.fits', '.h5')
    rrdesi_argument_zbest = output_argument.replace('coadd', 'zbest')

    print('srun -N 1 -n 32 -c 2 rrdesi_mpi -o {} -z {} {}'.format(rrdesi_argument_redrock, rrdesi_argument_zbest, output_argument))

