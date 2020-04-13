from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

n_plots = 100

for band in ['g', 'r', 'z']:

    if band=='g' or band=='r':
        ccd_dr9 = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr9sv/survey-ccds-90prime-dr9-cut.fits.gz')
    else:
        ccd_dr9 = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr9sv/survey-ccds-mosaic-dr9-cut.fits.gz')
    print(len(ccd_dr9))
    mask = ccd_dr9['ccd_cuts']==0
    ccd_dr9 = ccd_dr9[mask]
    print(len(ccd_dr9))

    ramin, ramax, decmin, decmax = 216, 230, 41, 46
    mask = (ccd_dr9['ra']>ramin) & (ccd_dr9['ra']<ramax) & (ccd_dr9['dec']>decmin) & (ccd_dr9['dec']<decmax)
    ccd_dr9 = ccd_dr9[mask]
    print(len(ccd_dr9))

    mask = ccd_dr9['filter']==band
    ccd_dr9_1 = ccd_dr9[mask]
    print(len(ccd_dr9_1))

    # random downsampling
    np.random.seed(123)
    idx = np.sort(np.random.choice(len(ccd_dr9_1), size=n_plots, replace=False))
    ccd_dr9_1 = ccd_dr9_1[idx]

    image_filename_dr8 = []
    if band=='g' or band=='r':
        ccd_dr8 = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/survey-ccds-90prime-dr8.fits.gz')
    else:
        ccd_dr8 = Table.read('/global/project/projectdirs/cosmo/data/legacysurvey/dr8/survey-ccds-mosaic-dr8.fits.gz')
    mask = np.in1d(ccd_dr8['expnum'], ccd_dr9_1['expnum'])
    ccd_dr8 = ccd_dr8[mask]
    for expnum in ccd_dr9_1['expnum']:
        mask = ccd_dr8['expnum']==expnum
        image_filename_dr8.append(ccd_dr8['image_filename'][mask][0].strip())

    ccd_dr9_1['median_dr9_'+band] = -99
    ccd_dr9_1['nmad_dr9_'+band] = -99
    ccd_dr9_1['median_dr8_'+band] = -99
    ccd_dr9_1['nmad_dr8_'+band] = -99

    for index in range(len(ccd_dr9_1)):

        ############### DR9 ###############

        photom_path = os.path.join('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/zpts', ccd_dr9_1['image_filename'][index].strip().replace('.fits.fz', '-photom.fits'))
        annotated_path = os.path.join('/global/cfs/cdirs/cosmo/work/legacysurvey/dr9/zpts', ccd_dr9_1['image_filename'][index].strip().replace('.fits.fz', '-annotated.fits'))

        photom = Table.read(photom_path)
        annotated = Table.read(annotated_path)

        if not np.all(annotated['zpt']==annotated['zpt'][0]):
            raise ValueError
        zpt = annotated['zpt'][0]

        mask = (photom['legacy_survey_mag']>16.5) & (photom['legacy_survey_mag']<18.5)
        mag_diff_median = np.median(photom['instpsfmag'][mask]+zpt-photom['legacy_survey_mag'][mask])
        mag_diff_nmad = nmad(photom['instpsfmag'][mask]+zpt-photom['legacy_survey_mag'][mask])
        ccd_dr9_1['median_dr9_'+band][index] = mag_diff_median
        ccd_dr9_1['nmad_dr9_'+band][index] = mag_diff_nmad

        # plt.figure(figsize=(6, 6))
        # plt.plot(photom['legacy_survey_mag'][~mask], photom['instpsfmag'][~mask]+zpt-photom['legacy_survey_mag'][~mask], '.', ms=1)
        # plt.plot(photom['legacy_survey_mag'][mask], photom['instpsfmag'][mask]+zpt-photom['legacy_survey_mag'][mask], '.', ms=1)
        # plt.xlabel('legacy_survey_mag')
        # plt.ylabel('(instpsfmag + zpt) - legacy_survey_mag')
        # plt.title('DR9 '+os.path.basename(photom_path).replace('-photom.fits', '')+'\n'+'Median = {:.4f}'.format(mag_diff_median)+',  $\sigma_{NMAD}$'+' = {:.4f}'.format(mag_diff_nmad))
        # plt.axis([14, 22, -0.2, 0.2])
        # plt.grid()
        # plt.savefig('/global/project/projectdirs/desi/www/users/rongpu/plots/dr9dev/ps1/mag_diff-{}-{}-dr9.png'.format(band, index))
        # plt.close()
        # # plt.show()

        ############### DR8 ###############

        download_dir = '/global/cscratch1/sd/rongpu/temp/dr8_photom'
        photom_path = os.path.join(download_dir, os.path.basename(image_filename_dr8[index].replace('.fits.fz', '-star-photom.fits')))
        # annotated_path = os.path.join(download_dir, os.path.basename(image_filename_dr8[index].replace('.fits.fz', '-star-annotated.fits')))

        if not os.path.exists(os.path.dirname(photom_path)):
            os.makedirs(os.path.dirname(photom_path))
        if (not os.path.isfile(photom_path)) or (os.stat(photom_path).st_size==0):
            if band=='g' or band=='r':
                url = 'https://faun.rc.fas.harvard.edu/eschlafly/dr8_zpt/bok/'+os.path.basename(photom_path)
            else:
                url = 'https://faun.rc.fas.harvard.edu/eschlafly/dr8_zpt/mosaic/'+os.path.basename(photom_path)
            cmd = 'wget -O '+photom_path+' \"'+url+'\"'
            print(cmd)
            download_status = os.system(cmd)
            if download_status!=0:
                continue
            # url = 'https://faun.rc.fas.harvard.edu/eschlafly/dr8_zpt/bok/'+os.path.basename(annotated_path)
            # cmd = 'wget -O '+annotated_path+' \"'+url+'\"'
            # print(cmd)
            # os.system(cmd)

        photom = Table.read(photom_path)
        # annotated = Table.read(annotated_path)

        # if not np.all(annotated['zpt']==annotated['zpt'][0]):
        #     raise ValueError
        # zpt = annotated['zpt'][0]

        mask = (photom['legacy_survey_mag']>16.5) & (photom['legacy_survey_mag']<18.5)
        mag_diff_median = np.median(photom['instpsfmag'][mask]+zpt-photom['legacy_survey_mag'][mask])
        mag_diff_nmad = nmad(photom['instpsfmag'][mask]+zpt-photom['legacy_survey_mag'][mask])
        ccd_dr9_1['median_dr8_'+band][index] = mag_diff_median
        ccd_dr9_1['nmad_dr8_'+band][index] = mag_diff_nmad

        # plt.figure(figsize=(6, 6))
        # plt.plot(photom['legacy_survey_mag'][~mask], photom['instpsfmag'][~mask]+zpt-photom['legacy_survey_mag'][~mask], '.', ms=1)
        # plt.plot(photom['legacy_survey_mag'][mask], photom['instpsfmag'][mask]+zpt-photom['legacy_survey_mag'][mask], '.', ms=1)
        # plt.xlabel('legacy_survey_mag')
        # plt.ylabel('(instpsfmag + zpt) - legacy_survey_mag')
        # plt.title('DR8 '+os.path.basename(photom_path).replace('-star-photom.fits', '')+'\n'+'Median = {:.4f}'.format(mag_diff_median)+',  $\sigma_{NMAD}$'+' = {:.4f}'.format(mag_diff_nmad))
        # plt.axis([14, 22, -0.2, 0.2])
        # plt.grid()
        # plt.savefig('/global/project/projectdirs/desi/www/users/rongpu/plots/dr9dev/ps1/mag_diff-{}-{}-dr8.png'.format(band, index))
        # plt.close()
        # # plt.show()
