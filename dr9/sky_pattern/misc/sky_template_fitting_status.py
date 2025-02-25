from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

image_dir = '/global/project/projectdirs/cosmo/staging'
blob_dir = '/global/cfs/cdirs/desi/users/rongpu/dr9/decam_ccd_blob_mask'
surveyccd_path = '/global/project/projectdirs/cosmo/work/legacysurvey/dr9/survey-ccds-decam-dr9.fits.gz'
template_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_templates_v2/'
skyscale_dir = '/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/sky_scales_v2/'

skyrun = Table.read('/global/cscratch1/sd/rongpu/temp/skyrunsgoodcountexpnumv48dr8.fits')
print(len(skyrun))

sky_path_list = glob.glob(os.path.join(template_dir, '*.fits.fz'))
print(len(sky_path_list))

# The file should be at least 5 hours old to ensure it's not being written
for sky_path in sky_path_list:
    time_modified = os.path.getmtime(sky_path)
    if (time.time() - time_modified)/3600 < 5:
        sky_path_list.remove(sky_path)

# #################################### Exclude z band ####################################
# sky_path_list = sky_path_list[::-1]
# mask = np.array(['_z_' in sky_path for sky_path in sky_path_list])
# sky_path_list = np.array(sky_path_list)[~mask]
# ########################################################################################

run_list = [int(fn[len(os.path.join(template_dir, 'sky_templates_'))+1:-8]) for fn in sky_path_list]
total_exposures = np.sum(np.in1d(skyrun['run'], run_list))

total_exposures_done = 0

for sky_path in sky_path_list:
    
    run = int(sky_path[len(os.path.join(template_dir, 'sky_templates_'))+1:-8])

    # Get run info
    mask = skyrun['run']==run
    n_exposure = np.sum(mask)
    band = skyrun['filter'][mask][0]
    mjd_span = skyrun['mjd_obs'][mask].max() - skyrun['mjd_obs'][mask].min()
    
    skyrun_idx = np.where(mask)[0]
    # print('\nrun {}, {} exposures'.format(run, len(skyrun_idx)))
    
    # Loop over exposures in the run
    for jj, index in enumerate(skyrun_idx):

        image_filename = skyrun['image_filename'][index].strip()
        expnum = skyrun['expnum'][index]

        skyscale_path = os.path.join(skyscale_dir, image_filename.replace('.fits.fz', '-skyscale.txt'))
        if os.path.isfile(skyscale_path):
            total_exposures_done += 1

print(total_exposures_done, 'exposures done')
print(total_exposures - total_exposures_done, 'exposures left')
print(total_exposures, 'total exposures')

############################### Save the per-run status in a table ###############################

run_status = Table()
run_status['run'] = np.unique(skyrun['run'])
run_status['filter'] = ' '
run_status['done'] = False
print(len(run_status['run']))

for run in run_list:
    
    # Get run info
    mask = skyrun['run']==run
    band = skyrun['filter'][mask][0]
    
    run_status_index = np.where(run_status['run']==run)
    run_status['filter'][run_status_index] = band

    skyrun_idx = np.where(mask)[0]
    # print('\nrun {}, {} exposures'.format(run, len(skyrun_idx)))
    
    # Loop over exposures in the run
    exposures_done = 0
    for jj, index in enumerate(skyrun_idx):

        image_filename = skyrun['image_filename'][index].strip()
        expnum = skyrun['expnum'][index]

        skyscale_path = os.path.join(skyscale_dir, image_filename.replace('.fits.fz', '-skyscale.txt'))
        if os.path.isfile(skyscale_path):
            exposures_done += 1
    
    if exposures_done==len(skyrun_idx):
        run_status['done'][run_status_index] = True
        
run_status.write('/global/cscratch1/sd/rongpu/dr9dev/sky_pattern/fitting_status.fits', overwrite=True)

print('Templates done:', np.sum(run_status['filter']!=' '))
print('Fitting done:', np.sum(run_status['done']))
for band in ['g', 'r', 'z']:
    mask = (run_status['filter']==band)
    mask1 = mask & (run_status['done'])
    print(np.sum(mask), np.sum(mask1))