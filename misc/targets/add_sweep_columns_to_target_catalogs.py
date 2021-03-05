# Add sweep, photo-z and stellar mass columns
# Example:
# python add_sweep_columns_to_target_catalogs.py sv1 LRG south

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio

from multiprocessing import Pool

from desitarget.targets import decode_targetid, encode_targetid

# Snippets taken from desitarget

def decode_sweep_name(sweepname):
    sweepname = os.path.basename(sweepname)

    ramin, ramax = float(sweepname[6:9]), float(sweepname[14:17])
    decmin, decmax = float(sweepname[10:13]), float(sweepname[18:21])

    if sweepname[9] == 'm':
        decmin *= -1
    if sweepname[17] == 'm':
        decmax *= -1

    return [ramin, ramax, decmin, decmax]


def is_in_box(objs, radecbox, ra_col='RA', dec_col='DEC'):

    ramin, ramax, decmin, decmax = radecbox

    # ADM check for some common mistakes.
    if decmin < -90. or decmax > 90. or decmax <= decmin or ramax <= ramin:
        msg = "Strange input: [ramin, ramax, decmin, decmax] = {}".format(radecbox)
        raise ValueError(msg)

    ii = ((objs[ra_col] >= ramin) & (objs[ra_col] < ramax)
          & (objs[dec_col] >= decmin) & (objs[dec_col] < decmax))

    return ii


n_processess = 32

ls_columns = ['FITBITS']

data_dir = '/global/cscratch1/sd/rongpu/target/catalogs/dr9.0/0.49.0'
stellar_mass_dir = '/global/cfs/cdirs/desi/users/rongpu/ls_dr9.0_photoz/stellar_mass'

# program: "main" or "sv1"
# target_class: "LRG", "ELG", "QSO" or "BGS_ANY"
# field: "north" or "south"
program, target_class, field = str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3])
program = program.lower()
target_class = target_class.upper()
field = field.lower()

print(program, target_class, field)

cat_basic_path = os.path.join(data_dir, 'dr9_{}_{}_{}_0.49.0_basic.fits'.format(program, target_class.lower(), field))
more_path = os.path.join(data_dir, 'dr9_{}_{}_{}_0.49.0_more_1.fits'.format(program, target_class.lower(), field))

if os.path.isfile(more_path):
    sys.exit('File already exist: '+more_path)

cat_basic = Table(fitsio.read(cat_basic_path, columns=['RA', 'DEC', 'TARGETID']))

# #########################################################################################
# cat_basic = cat_basic[:len(cat_basic)//50]
# #########################################################################################

sweep_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/{}/sweep/9.0'.format(field)
sweep_fn_list = np.array(sorted(glob.glob(os.path.join(sweep_dir, '*.fits'))))

sweep_radec_list = [decode_sweep_name(sweep_fn) for sweep_fn in sweep_fn_list]
mask = np.array([np.any(is_in_box(cat_basic, sweep_radec)) for sweep_radec in sweep_radec_list])
print(np.sum(mask), len(mask))
sweep_fn_list = sweep_fn_list[mask]


def get_sweep_columns(sweep_fn):

    cat = Table(fitsio.read(sweep_fn, columns=['OBJID', 'BRICKID', 'RELEASE']))
    targetid = encode_targetid(cat['OBJID'], cat['BRICKID'], cat['RELEASE'])
    idx = np.where(np.in1d(targetid, cat_basic['TARGETID']))[0]
    if len(idx)==0:
        return None
    targetid = targetid[idx]
    cat = Table(fitsio.read(sweep_fn, rows=idx, columns=ls_columns))
    cat['TARGETID'] = targetid
    pz_fn = sweep_fn.replace('sweep/9.0/', 'sweep/9.0-photo-z/').replace('.fits', '-pz.fits')
    pz = Table(fitsio.read(pz_fn, rows=idx))
    pz.remove_columns(['OBJID', 'BRICKID', 'RELEASE'])
    cat = hstack([cat, pz], join_type='exact')

    # Add stellar mass
    stellar_mass_path = os.path.join(stellar_mass_dir, field, os.path.basename(sweep_fn).replace('.fits', '_stellar_mass.npy'))
    cat['stellar_mass'] = np.load(stellar_mass_path)[idx]

    return cat


if __name__ == '__main__':

    print('Start!')
    time_start = time.time()

    # start multiple worker processes
    with Pool(processes=n_processess) as pool:
        res = pool.map(get_sweep_columns, np.unique(sweep_fn_list))

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

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
