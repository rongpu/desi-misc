# Merge spectra of a list of TARGETIDs from uniqpix-based coadd files into a single coadd file.

import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp
import multiprocessing
from multiprocessing import Pool

from desiutil.log import get_logger
from desispec.io import read_spectra, write_spectra, findfile
from desispec.spectra import stack

log = get_logger()

#- Python 3.14 defaults to forkserver, which cannot find functions defined in an interactive __main__
multiprocessing.set_start_method('fork', force=True)

#----------------------------------------------------------------------------

fiberassign_file = '/global/cfs/cdirs/desi/survey/fiberassign/special/tertiary/0026/tertiary-targets-0026-assign.fits'  # path to a file with the TARGETIDs; see load_targetids()

specprod = 'matterhorn'
survey = 'special'
program = 'dark'

zpix_file = '/global/cfs/cdirs/desi/spectro/redux/matterhorn/zcatalog/v2/special/zpix-special-dark.fits'  # path to zpix catalog with TARGETID and UNIQPIX

outfile = '/pscratch/sd/r/rongpu/desi2/uniqpix_coadds/coadd-tertiary26-laelbg.fits'

parallel = True
nproc = 128

#----------------------------------------------------------------------------


def read_one_file(args):
    """Read the requested targets from one uniqpix coadd file.

    Args:
        args: tuple (uniqpix, survey, program, targetids)

    Returns:
        Spectra object for the targets found in the file, or None if none found
    """
    uniqpix, survey, program, targetids = args
    specfile = findfile('coadd_pix', uniqpix=uniqpix, survey=survey, faprogram=program,
                        specprod=specprod, readonly=True)
    log.info(f'Reading {len(targetids)} targets from {specfile}')
    sp = read_spectra(specfile, targetids=targetids)

    if sp.num_targets() == 0:
        log.warning(f'No matching targets found in {specfile}')
        return None

    #- keep provenance since the merged file mixes pixels, surveys, and programs
    sp.fibermap['UNIQPIX'] = uniqpix
    sp.fibermap['SURVEY'] = survey
    sp.fibermap['PROGRAM'] = program

    return sp


fa = Table(fitsio.read(fiberassign_file))
zpix = Table(fitsio.read(zpix_file, columns=['TARGETID', 'UNIQPIX']))
zpix['SURVEY'] = survey
zpix['PROGRAM'] = program
print(len(fa))
print(len(zpix))
print()

# Sanity checks
mask0 = np.isin(zpix['TARGETID'], fa['TARGETID'])
print(np.sum(mask0))
mask = fa['NASSIGN']>0
fa = fa[mask]
print(len(fa))
mask1 = np.isin(zpix['TARGETID'], fa['TARGETID'])
print(np.sum(mask1))
assert np.sum(mask0)==np.sum(mask1)
assert np.sum(mask0)==len(fa)
print()

# Remove "ELG" spectra
mask = fa['TERTIARY_TARGET']!='ELG'
fa = fa[mask]
print(len(fa))

#- map TARGETID -> (UNIQPIX, SURVEY, PROGRAM); a TARGETID can have several rows
zpix = zpix[np.isin(zpix['TARGETID'], fa['TARGETID'])]
print(len(zpix))

missing = np.setdiff1d(fa['TARGETID'], zpix['TARGETID'])
if len(missing) > 0:
    log.warning(f'{len(missing)} TARGETIDs not found in {zpix_file}')

#- one task per coadd file, so each file is only read once
arglist = []
for zz in zpix.group_by(['UNIQPIX', 'SURVEY', 'PROGRAM']).groups:
    arglist.append((int(zz['UNIQPIX'][0]), str(zz['SURVEY'][0]), str(zz['PROGRAM'][0]),
                    np.array(zz['TARGETID'])))
log.info(f'{len(zpix)} spectra in {len(arglist)} coadd files')

if parallel:
    with Pool(min(nproc, len(arglist))) as pool:
        spectra = pool.map(read_one_file, arglist)
else:
    spectra = [read_one_file(args) for args in arglist]
spectra = [sp for sp in spectra if sp is not None]

if len(spectra) == 0:
    raise RuntimeError('No matching targets found in any file')

merged = stack(spectra)
write_spectra(outfile, merged)
log.info(f'Wrote {merged.num_spectra()} spectra to {outfile}')

