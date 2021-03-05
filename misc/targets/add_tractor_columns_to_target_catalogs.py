# Add tractor columns
# Example:
# python add_tractor_columns_to_target_catalogs.py sv1 LRG south

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack
import fitsio

from multiprocessing import Pool

from desitarget.targets import decode_targetid, encode_targetid

n_processess = 32

data_dir = '/global/cscratch1/sd/rongpu/target/catalogs/dr9.0/0.49.0'

# Add more columns from tractor catalogs
tractor_columns = ['fitbits', 'nea_g', 'nea_r', 'nea_z', 'blob_nea_g', 'blob_nea_r', 'blob_nea_z']

# program: "main" or "sv1"
# target_class: "LRG", "ELG", "QSO" or "BGS_ANY"
# field: "north" or "south"
program, target_class, field = str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3])
program = program.lower()
target_class = target_class.upper()
field = field.lower()

if field=='south':
    photsys = 'S'
elif field=='north':
    photsys = 'N'

print(program, target_class, field)

cat_basic_path = os.path.join(data_dir, 'dr9_{}_{}_{}_0.49.0_basic.fits'.format(program, target_class.lower(), field))
more_path = os.path.join(data_dir, 'dr9_{}_{}_{}_0.49.0_more_0.fits'.format(program, target_class.lower(), field))

if os.path.isfile(more_path):
    sys.exit('File already exist: '+more_path)

cat_basic = Table(fitsio.read(cat_basic_path, columns=['TARGETID']))
objid_all, brickid_all, release_all = decode_targetid(cat_basic['TARGETID'])[:3]

bricks = Table.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/survey-bricks.fits.gz')

tractor_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/tractor'.format(field)


def get_tractor_columns(brickid):

    brickname = bricks['BRICKNAME'][bricks['BRICKID']==brickid][0]
    tractor_path = os.path.join(tractor_dir, brickname[:3], 'tractor-'+brickname+'.fits')
    tractor_objid = fitsio.read(tractor_path, columns=['objid'])['objid']
    mask = brickid_all==brickid
    idx = np.where(np.in1d(tractor_objid, objid_all[mask]))[0]
    cat = Table(fitsio.read(tractor_path, columns=tractor_columns+['objid', 'brickid', 'release'], rows=idx))
    cat['TARGETID'] = encode_targetid(cat['objid'], cat['brickid'], cat['release'])
    cat.remove_columns(['objid', 'brickid', 'release'])

    if len(cat)!=np.sum(brickid_all==brickid):
        print(len(cat), np.sum(brickid_all==brickid))
        raise ValueError('different catalog length')

    return cat


if __name__ == '__main__':

    print('Start!')
    time_start = time.time()

    # start multiple worker processes
    with Pool(processes=n_processess) as pool:
        res = pool.map(get_tractor_columns, np.unique(brickid_all))

    # # Remove None elements from the list
    # for index in range(len(res)-1, -1, -1):
    #     if res[index] is None:
    #         res.pop(index)

    cat_more = vstack(res, join_type='exact')
    if len(cat_more)!=len(cat_basic):
        print(len(cat_more), len(cat_basic))
        raise ValueError('different catalog length')

    # Here matching cat_more to cat_basic
    t1_reverse_sort = np.array(cat_basic['TARGETID']).argsort().argsort()
    cat_more = cat_more[np.argsort(cat_more['TARGETID'])[t1_reverse_sort]]
    if not np.all(cat_more['TARGETID']==cat_basic['TARGETID']):
        raise ValueError('different targetid')
    cat_more.remove_column('TARGETID')

    cat_more.write(more_path)

    print(time.strftime("%H:%M:%S", time.gmtime(time.time() - time_start)))
