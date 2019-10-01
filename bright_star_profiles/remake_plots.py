# Run on NERSC

from __future__ import division, print_function
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
import sys, os, glob, time, warnings, gc
import healpy as hp

params = {'legend.fontsize': 'large',
         'axes.labelsize': 'large',
         'axes.titlesize':'large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large',
         'figure.facecolor':'w'} 
plt.rcParams.update(params)

###############################################################

download_dir = '/global/cscratch1/sd/rongpu/temp/dr8_cutouts'

# # DECaLS NGC
# ramin, ramax, decmin, decmax = 160, 210, -5, 20
# field = 'south'
# output_name = 'decals_ngc'

# # DECaLS SGC
# ramin, ramax, decmin, decmax = 350, 35, 6, 25
# field = 'south'
# output_name = 'decals_sgc'

# # DES
# ramin, ramax, decmin, decmax = 24, 44, -25, 5
# ramin2, ramax2, decmin2, decmax2 = 1, 44, -6, 5
# field = 'south'
# output_name = 'des'

# BASS/MZLS
ramin, ramax, decmin, decmax = 150, 220, 35, 65
field = 'north'
output_name = 'bassmzls'

print(field, output_name)

###############################################################

# img_type = '-resid'
img_type = ''
size_str = '512'

# narrow bins of LS magnitude
ls_mag_bins = [10.5, 11.75, 13.0, 14.25, 15.5]

###################################################################################################

for band in ['g', 'r', 'z']:

    gaia_output = Table.read('data/individual_profiles_{}_{}_{}.fits'.format(field, output_name, band))
    profile = Table.read('data/average_profiles_{}_{}_{}.txt'.format(field, output_name, band), format='ascii.commented_header')

    for index in range(len(ls_mag_bins)):
        mask = gaia_output['ls_'+band+'_bin']==ls_mag_bins[index]
        plt.figure(figsize=(8, 6))
        plt.loglog(gaia_output['radius'][mask].T, gaia_output['flux'][mask].T, lw=1, alpha=0.1, c='C0')
        plt.loglog(profile['radius_{}_{:.2f}'.format(band, ls_mag_bins[index])], profile['flux_{}_{:.2f}'.format(band, ls_mag_bins[index])], lw=4, alpha=0.5, c='C1')
        # plt.errorbar(profile['radius_{}_{:.2f}'.format(band, ls_mag_bins[index])], profile['flux_{}_{:.2f}'.format(band, ls_mag_bins[index])], yerr=flux_spread_in_bin[index],
        #              lw=1, alpha=1, c='C1', zorder=5)
        plt.axis([.5, 70, .001, 200])
        plt.xlabel('Radius (arcsec)')
        plt.ylabel('SB (nmgy/sq.arcsec.)')
        plt.grid(alpha=0.5)
        plt.title('{} {} {}mag = {:.2f}'.format(field, output_name, band, ls_mag_bins[index]))
        plt.savefig('plots/{}_{}_{}mag_{:.2f}.png'.format(field, output_name, band, ls_mag_bins[index]))
        plt.close()

    plt.figure(figsize=(11, 8))
    for index in range(len(ls_mag_bins)):
        # normalize the flux to 13th magnitude stars
        norm = 10**((ls_mag_bins[index]-13)/2.5)
        plt.loglog(profile['radius_{}_{:.2f}'.format(band, ls_mag_bins[index])], profile['flux_{}_{:.2f}'.format(band, ls_mag_bins[index])]*norm, lw=1.5, alpha=1., 
                   label='{}mag = {:.2f}'.format(band, ls_mag_bins[index]), c='C'+str(index))
        # plt.errorbar(profile['radius_{}_{:.2f}'.format(band, ls_mag_bins[index])], profile['flux_{}_{:.2f}'.format(band, ls_mag_bins[index])]*norm, yerr=flux_spread_in_bin[index],
        #              lw=1, alpha=1, c='C'+str(index), zorder=5)
    plt.title('{} {} {}-band'.format(field, output_name, band))
    plt.axis([.5, 100, .0002, 1000])
    plt.grid(alpha=0.5)
    plt.xlabel('Radius (arcsec)')
    plt.ylabel('SB (a.u.)')
    plt.legend()
    plt.savefig('plots/{}_{}_{}mag_average.png'.format(field, output_name, band))
    plt.close()

    plt.figure(figsize=(11, 8))
    for index in range(len(ls_mag_bins)):
        # normalize the flux to 13th magnitude stars
        norm = 10**((ls_mag_bins[index]-13)/2.5)
        plt.loglog(profile['radius_{}_{:.2f}'.format(band, ls_mag_bins[index])], profile['flux_{}_{:.2f}'.format(band, ls_mag_bins[index])]*norm, lw=1.5, alpha=1., 
                   label='{}mag = {:.2f}'.format(band, ls_mag_bins[index]), c='C'+str(index))
    plt.title('{} {} {}-band'.format(field, output_name, band))
    plt.axis([2, 20, 0.01, 20])
    plt.grid(alpha=0.5)
    plt.xlabel('Radius (arcsec)')
    plt.ylabel('SB (a.u.)')
    plt.legend()
    plt.savefig('plots/{}_{}_{}mag_average_zoomin.png'.format(field, output_name, band))
    plt.close()

    ###################################################################################################

