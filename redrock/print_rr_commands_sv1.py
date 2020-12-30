# Get the a list of science exposures for a specficic date and tile
# Print desi_coadd_spectra commands for coadds of single or multiple exposures
# Print redrock commands

# master branch:
# source /global/cfs/cdirs/desi/software/desi_environment.sh master
# export REDUXDIR=/global/cfs/cdirs/desi/spectro/redux/daily
# export OUTDIR=$SCRATCH/desi/sv1/deep_coadd

from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

redux_dir = '/global/cfs/cdirs/desi/spectro/redux/daily'
output_dir = '/global/cscratch1/sd/rongpu/desi/sv1/deep_coadd'
# redux_dir = os.getenv('REDUXDIR')
# output_dir = os.getenv('OUTDIR')

obsdate_list = ['20201214', '20201215', '20201216', '20201217', '20201218', '20201219', '20201220', '20201221']
n_exp = np.inf        # number of exposures in a coadded; 1 for single-exposure coadd;
                      # use np.inf to coadd all exposures

overwrite = False

# tileid_list = None  # no restriction on tiles
tileid_list = [80605, 80607, 80609, 80620, 80622]

################################## Get list of exposures ##################################

exposure_dir_list = []
for obsdate in obsdate_list:
    exposure_dir_list += glob.glob(os.path.join(redux_dir, 'exposures', obsdate, '*'))

cframe_list = []

# Get a list of all science exposures
for exposure_dir in exposure_dir_list:

    cframe_list_tmp = glob.glob(os.path.join(exposure_dir, 'cframe-*'))

    if len(cframe_list_tmp)>0:

        if tileid_list is None:
            cframe_list += cframe_list_tmp
        else:
            # only need to check one cframe file in the exposure
            with fitsio.FITS(cframe_list_tmp[0]) as f:
                if f[0].read_header()['TILEID'] in tileid_list:
                    cframe_list += cframe_list_tmp

cframe_list = sorted(cframe_list)

# Gather exposure/petal information
cframes = Table()
cframes['cframe'] = np.array(cframe_list)
cframes['tileid'] = np.zeros(len(cframes), dtype=int)
cframes['expid'] = np.zeros(len(cframes), dtype=int)
cframes['camera'] = ' '
cframes['petal_loc'] = -1 * np.ones(len(cframes), dtype=np.int32)
for index, cframe in enumerate(cframes['cframe']):
    # cframes['camera'][index] = cframes[cframe.find('/cframe-')+len('/cframe-'):][:2][0]
    # cframes['petal_loc'][index] = cframes[cframe.find('/cframe-')+len('/cframe-'):][:2][1]
    with fitsio.FITS(cframe) as f:
        header = f[0].read_header()
        cframes['tileid'][index] = header['TILEID']
        cframes['expid'][index] = header['EXPID']
        cframes['camera'][index] = header['CAMERA'].strip()[0]
        cframes['petal_loc'][index] = int(header['CAMERA'].strip()[1])

# Sanity check: each petal must have three cframe files
for expid in np.unique(cframes['expid']):
    mask_expid = cframes['expid']==expid
    for petal_loc in range(10):
        mask = mask_expid & (cframes['petal_loc']==petal_loc)
        if (np.sum(mask)>0) & (np.sum(mask)!=3):
            raise ValueError('EXPID {} PETAL_LOC {} has only {} cframes files'.format(expid, petal_loc, np.sum(mask)))

################################## Print desi_coadd_spectra commands ##################################

print('\n################# desi_coadd_spectra commands: #################\n')

output_argument_list = []

for tileid in np.unique(cframes['tileid']):

    for petal_loc in range(10):

        mask = (cframes['tileid']==tileid) & (cframes['petal_loc']==petal_loc)
        mask &= (cframes['camera']=='b') # choose one camera for simplicity

        if (np.sum(mask)<n_exp) and (n_exp!=np.inf):
            # print('\n# Not enough exposures in TILEID {} PETAL_LOD {} for one coadd\n'.format(tileid, petal_loc))
            continue

        cframe1 = cframes[mask]

        # Skip the exposures that are that do not make the split
        if n_exp!=np.inf:
            cframe1 = cframe1[:len(cframe1)-len(cframe1)%(n_exp)]
            nsplit = len(cframe1)//(n_exp)
            subset_split = np.split(np.arange(len(cframe1)), nsplit)
        else:
            subset_split = [np.arange(len(cframe1))]

        for subset_index in range(len(subset_split)):

            subset = cframe1[subset_split[subset_index]]
            input_argument = ''

            for index in range(len(subset)):
                exposure_dir = os.path.dirname(subset['cframe'][index])
                input_argument += os.path.join(exposure_dir.replace(redux_dir, '$REDUXDIR'), 'cframe-[brz]{}-*.fits ').format(petal_loc)

            # print(input_argument)

            if n_exp==np.inf:
                if len(obsdate_list)==1:
                    # same filename format as andes
                    output_argument = os.path.join('$OUTDIR', str(tileid), 'coadd-{}-{}-{}.fits'.format(petal_loc, tileid, obsdate_list[0]))
                else:
                    output_argument = os.path.join('$OUTDIR', str(tileid), 'coadd-{}-{}.fits'.format(petal_loc, tileid))
            elif n_exp==1:
                exposure = os.path.basename(exposure_dir)
                output_argument = os.path.join('$OUTDIR', str(tileid), 'coadd-{}-{}.fits'.format(petal_loc, exposure))
            else:
                output_argument = os.path.join('$OUTDIR', str(tileid), 'coadd-{}-{}exp-subset-{}.fits'.
                format(petal_loc, n_exp, subset_index))

            if os.path.isfile(os.path.expandvars(output_argument)) and (not overwrite):
                print('\nWarninig: {} already exists!\n'.format(output_argument))
                continue

            output_argument_list.append(output_argument)

            print('time desi_coadd_spectra --coadd-cameras -i {} -o {}'.format(input_argument, output_argument))
            # print('time desi_coadd_spectra --inframes {} --outfile {}'.format(input_argument, output_argument))

################################## Print redrock commands ##################################

print('\n################# redrock commands: #################\n')

for output_argument in output_argument_list:

    rrdesi_argument_redrock = output_argument.replace('coadd', 'redrock')
    rrdesi_argument_zbest = output_argument.replace('coadd', 'zbest')

    print('srun -N 1 -n 32 -c 2 rrdesi_mpi {} -o {} -z {}'.format(output_argument, rrdesi_argument_redrock, rrdesi_argument_zbest))

