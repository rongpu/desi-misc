# Ignore z camera for redshift fitting by setting Z_IVAR to 0
# ln -s /pscratch/sd/r/rongpu/desi2/tertiary49/ignore_z_camera/* /global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/ignore_z_camera/

import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

# fns = glob.glob('/global/cfs/cdirs/desi/spectro/redux/darknight/tiles/cumulative/83577/20250430/coadd-*.fits')
# fns.sort()

# output_dir = '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary47/remove_z_camera/'

tile_list = [83579, 83580, 83581, 83582]
lastnight_list = [20260116, 20260115, 20260116, 20260116]
fns = []
for tileid, lastnight in zip(tile_list, lastnight_list):
    print(tileid, lastnight)
    fns += glob.glob('/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/{}/{}/coadd-*thru{}.fits'.format(tileid, lastnight, lastnight))
print(len(fns))

output_dir = '/pscratch/sd/r/rongpu/desi2/tertiary49/ignore_z_camera/'

for fn in fns:
    fn_output = os.path.join(output_dir, os.path.basename(fn))
    fits = fitsio.FITS(fn)
    with fitsio.FITS(fn_output, 'rw') as f:
        for index in range(len(fits)):
            extname = fits[index].get_extname()
            data, header = fitsio.read(fn, ext=index, header=True)
            if extname=='Z_IVAR':
                data *= 0.
                assert np.all(data==0)
            f.write(data, extname=extname, header=header)
    fits.close()
