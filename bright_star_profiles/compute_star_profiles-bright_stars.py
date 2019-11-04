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

# DECaLS NGC
ramin, ramax, decmin, decmax = 160, 210, -5, 30
field = 'south'
output_name = 'decals_ngc'

# # DECaLS SGC
# ramin, ramax, decmin, decmax = 350, 35, 6, 25
# field = 'south'
# output_name = 'decals_sgc'

# # DES
# ramin, ramax, decmin, decmax = 24, 44, -25, 5
# ramin2, ramax2, decmin2, decmax2 = 1, 44, -6, 5
# field = 'south'
# output_name = 'des'

# # BASS/MZLS
# ramin, ramax, decmin, decmax = 150, 220, 35, 65
# field = 'north'
# output_name = 'bassmzls'

print(field, output_name)

###############################################################

# img_type = '-resid'
img_type = ''
size_str = '512'

# narrow bins of LS magnitude
ls_mag_bins = [6.5]
ls_mag_bin_width = 1.
nsamp = 200
# plot the image and profiles of individual objects
individual_plot_q = False
# plot_dir = '/global/homes/r/rongpu/notebooks/star_profiles/plots/11/'
verbose = False


def plot_cutout(img, pixscale, vmin=-1, vmax=1, unit='arcsec'):
    if unit=='arcsec':
        extent = 0.5*pixscale*img.shape[0]*np.array([-1, 1, -1, 1])
    elif unit=='arcmin':
        extent = 0.5*pixscale*img.shape[0]*np.array([-1, 1, -1, 1])/60.
    elif unit=='deg':
        extent = 0.5*pixscale*img.shape[0]*np.array([-1, 1, -1, 1])/3600.
    else:
        raise ValueError('unrecognized unit')
    fig, ax = plt.subplots(figsize=(8, 8))
    dens = ax.imshow(img, aspect='equal', 
               cmap='seismic', extent=extent, vmin=vmin, vmax=vmax)
    ax.axvline(0, ls='--', lw=0.5, color='grey')
    ax.axhline(0, ls='--', lw=0.5, color='grey')
    fig.colorbar(dens, fraction=0.046, pad=0.04)
    return ax

def binned_stats(x, y, bins):
    nmad = lambda x: 1.4826*np.median(np.abs(x-np.median(x)))
    bin_center, bin_median, bin_scatter = np.zeros((3, len(bins)-1))
    for index in range(len(bins)-1):
        mask = (x>bins[index]) & (x<bins[index+1])
        bin_center[index] = np.median(x[mask])
        if np.sum(mask)>0:
            bin_median[index] = np.median(y[mask])
            bin_scatter[index] = nmad(y[mask])
        else:
            bin_median[index], bin_scatter[index] = np.nan, np.nan
    return bin_center, bin_median, bin_scatter


gaia_dir = '/project/projectdirs/cosmo/data/gaia/dr2/healpix'

gaia_nside = 32
gaia_npix = hp.nside2npix(gaia_nside)
gaia_hp_ra, gaia_hp_dec = hp.pix2ang(gaia_nside, np.arange(gaia_npix), nest=True, lonlat=True)
if ramin<ramax:
    mask = (gaia_hp_ra>ramin) & (gaia_hp_ra<ramax) & (gaia_hp_dec>decmin) & (gaia_hp_dec<decmax)
else:
    mask = ((gaia_hp_ra>ramin) | (gaia_hp_ra<ramax)) & (gaia_hp_dec>decmin) & (gaia_hp_dec<decmax)
if output_name=='des':
    mask |= (gaia_hp_ra>ramin2) & (gaia_hp_ra<ramax2) & (gaia_hp_dec>decmin2) & (gaia_hp_dec<decmax2)
gaia_list = np.where(mask)[0]
gaia = []
for hp_idx in gaia_list:
    gaia_fn = (5-len(str(hp_idx)))*'0'+str(hp_idx)
    tmp = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), colnames=['PHOT_G_MEAN_MAG'])
    # select GAIA_G<17.5 objects
    idx = np.where(tmp['PHOT_G_MEAN_MAG']<17.5)[0]
    tmp = fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), rows=idx)
    tmp = Table(tmp)
    gaia.append(tmp)
gaia = vstack(gaia)
print(len(gaia))

# convert column names to lower case
for col in gaia.colnames:
    gaia.rename_column(col, col.lower())
    
if ramin<ramax:
    mask = (gaia['ra']>ramin) & (gaia['ra']<ramax) & (gaia['dec']>decmin) & (gaia['dec']<decmax)
else:
    mask = ((gaia['ra']>ramin) | (gaia['ra']<ramax)) & (gaia['dec']>decmin) & (gaia['dec']<decmax)
if output_name=='des':
    mask |= (gaia['ra']>ramin2) & (gaia['ra']<ramax2) & (gaia['dec']>decmin2) & (gaia['dec']<decmax2)
gaia = gaia[mask]
print(len(gaia))

if field=='north':
    mask = gaia['dec']>32.375
else:
    mask = gaia['dec']<=32.375    
gaia = gaia[mask]
print(len(gaia))
    
# Remove the globular clusters and bad regions
mask = (gaia['ra']>198.0) & (gaia['ra']<198.5) & (gaia['dec']>18.0) & (gaia['dec']<18.4)
mask |=  (gaia['ra']>199.0) & (gaia['ra']<199.2) & (gaia['dec']>17.6) & (gaia['dec']<17.8)
mask |=  (gaia['ra']>182.4) & (gaia['ra']<182.65) & (gaia['dec']>18.44) & (gaia['dec']<18.64)
mask |= (gaia['ra']>13.0) & (gaia['ra']<13.4) & (gaia['dec']>-26.7) & (gaia['dec']<-26.45)
mask |= (gaia['ra']>14.5) & (gaia['ra']<15.5) & (gaia['dec']>-34.0) & (gaia['dec']<-33.2)
gaia = gaia[~mask]
print(len(gaia))

gaia['bp_rp'] = gaia['phot_bp_mean_mag'] - gaia['phot_rp_mean_mag']
gaia['pm'] = np.sqrt(gaia['pmra']**2 + gaia['pmdec']**2)

# Remove duplicates
if len(np.unique(gaia['source_id']))<len(gaia):
    print('Duplicates exist!')
    gaia.sort('source_id')
    mask = gaia['source_id'][1:]==gaia['source_id'][:-1]
    mask = np.concatenate([[False], mask])
    gaia = gaia[~mask]
    
# Remove duplicated_source==True
mask = gaia['duplicated_source'].copy()
if np.sum(mask)>0:
    print('{} objects with duplicated_source==True'.format(np.sum(mask)))
    gaia = gaia[~mask]
    print(len(gaia))

# Remove objects with invalid GAIA color
mask = ~np.isfinite(gaia['bp_rp'])
mask |= (gaia['phot_bp_mean_mag']==0) | (gaia['phot_rp_mean_mag']==0)
if np.sum(mask)>0:
    print('{} objects with invalid bp_rp'.format(np.sum(mask)))
    gaia = gaia[~mask]
    print(len(gaia))

# Apply GAIA-LS transformation
from scipy.interpolate import interp1d
transform_interp = {}
for band in ['g', 'r', 'z']:
    if field=='north':
        tmp = Table.read('data/gaia_bassmzls_transform.txt', format='ascii.commented_header')
    else:
        tmp = Table.read('data/gaia_decals_transform.txt', format='ascii.commented_header')
    transform_interp[band] = interp1d(tmp['bp_rp'], tmp['ls_'+band], bounds_error=False, fill_value='extrapolate', kind='linear')
    gaia['ls_'+band] = gaia['phot_g_mean_mag'] + transform_interp[band](gaia['bp_rp'])

# plt.figure(figsize=(15, 10))
# if ramin<ramax:
#     ra_plot = gaia['ra'].copy()
# else:
#     ra_plot = (gaia['ra']+180)%360 - 180
# plt.plot(ra_plot[::1], gaia['dec'][::1], '.', ms=0.1, alpha=0.5)
# plt.grid(alpha=0.5)
# plt.savefig('plots/more/{}_{}_coverage.png'.format(field, output_name))
# plt.close()

###################################################################################################

for band in ['g', 'r', 'z']:

    # The native pixel size is 0.262 in both north and south in the cutouts
    # pixscale_str = '0.262'
    pixscale_str = '1.048'
    radius_binsize = 1.

    pixscale = float(pixscale_str)

    radius_in_bin, flux_in_bin, flux_scatter_in_bin = [], [], []
    gaia_output = []

    for mag_index in range(len(ls_mag_bins)):

        ls_mag_min, ls_mag_max = ls_mag_bins[mag_index]-ls_mag_bin_width/2, ls_mag_bins[mag_index]+ls_mag_bin_width/2
        print('{} < {} < {}'.format(ls_mag_min, band, ls_mag_max))

        idx = np.where((gaia['ls_'+band]>ls_mag_min) & (gaia['ls_'+band]<ls_mag_max))[0]
        print(len(idx))
        
        # cut on proper motion
        mask = gaia['pm'][idx]<100
        idx = idx[mask]
        print(len(idx))
        
        if len(idx)>nsamp:
            np.random.seed(1)
            idx = np.sort(np.random.choice(idx, size=nsamp, replace=False))
        else:
            print('!!!Only {} stars are used!!!'.format(len(idx)))
        
        # if not os.path.exists(plot_dir):
        #     os.makedirs(plot_dir)

        radius_arr, flux_arr = [], []

        for index in idx:

            if verbose:
                print('phot_g_mean_mag = ', gaia['phot_g_mean_mag'][index])

            ra, dec = gaia['ra'][index], gaia['dec'][index]

            file_path = os.path.join(download_dir, field, '{}_{}{}_{}_{}.fits'.format(gaia['source_id'][index], band, img_type, pixscale_str, size_str))
            if not os.path.exists(os.path.dirname(file_path)):
                os.makedirs(os.path.dirname(file_path))

            if not os.path.isfile(file_path):
                url = 'http://legacysurvey.org/viewer/cutout.fits?ra={:f}&dec={:f}&layer=dr8-{}{}&pixscale={}&bands={}&size={}'.format(ra, dec, field, img_type, pixscale_str, band, size_str)
                cmd = 'wget -O '+file_path+' \"'+url+'\"'
                if verbose:
                    print(cmd)
                os.system(cmd)

            img_raw = fitsio.read(file_path)
            # normalize to actual flux per pixel
            img = img_raw / (pixscale**2) * (pixscale/0.262)**2

            grid = pixscale * np.linspace(-0.5*(img.shape[0]-1), 0.5*(img.shape[0]-1), img.shape[0])
            xx, yy = np.meshgrid(grid, grid)
            radius = np.sqrt(xx**2 + yy**2).flatten()
            max_radius = xx.flatten().max() # maximum radius before hitting the edge of the image
            nbins = int(np.floor(max_radius/radius_binsize)+1)

            # Remove masked pixels
            mask = img.flatten()!=0
            radius = radius[mask]
            flux = img.flatten()[mask]

            bin_center, bin_median, bin_scatter = binned_stats(radius, flux, bins=np.linspace(0., np.floor(max_radius), nbins))
            # normalize to the magnitude bin center
            bin_median = bin_median * 10**((gaia['ls_'+band][index]-ls_mag_bins[mag_index])/2.5)
            radius_arr.append(bin_center)
            flux_arr.append(bin_median)
            
            # if individual_plot_q:
            #     vrange = 0.5
            #     ax = plot_cutout(img, pixscale, vmin=-vrange, vmax=vrange)
            #     ax.set_title('['+band+'-band]  GAIA_G={:.4f}'.format(gaia['phot_g_mean_mag'][index]))
            #     plt.savefig(os.path.join(plot_dir, os.path.basename(file_path)[:-5]+'_image.png'))
            #     plt.close()

            #     plt.figure(figsize=(8, 6))
            #     plt.plot(radius, flux, '.', ms=0.5)
            #     plt.plot(bin_center, bin_median, c='C1')
            #     plt.errorbar(bin_center, bin_median, yerr=bin_scatter, lw=1, alpha=0.6, c='C1')
            #     plt.axis([0, 70, -1, 5])
            #     plt.axhline(0, lw=1, color='r')
            #     plt.grid(alpha=0.5)
            #     plt.title('['+band+'-band]  GAIA_G={:.4f}'.format(gaia['phot_g_mean_mag'][index]))
            #     plt.savefig(os.path.join(plot_dir, os.path.basename(file_path)[:-5]+'_profile.png'))
            #     plt.close()

            #     plt.figure(figsize=(8, 6))
            #     plt.loglog(radius, flux, '.', ms=0.5)
            #     plt.plot(bin_center, bin_median, c='C1')
            #     plt.errorbar(bin_center, bin_median, yerr=bin_scatter, lw=1, alpha=0.6, c='C1')
            #     plt.axis([.5, 70, .02, 200])
            #     plt.grid(alpha=0.5)
            #     plt.title('['+band+'-band]  GAIA_G={:.4f}'.format(gaia['phot_g_mean_mag'][index]))
            #     plt.savefig(os.path.join(plot_dir, os.path.basename(file_path)[:-5]+'_profile_log.png'))
            #     plt.close()
        
        radius_arr = np.array(radius_arr)
        flux_arr = np.array(flux_arr)
        
        gaia_output_tmp = gaia[idx].copy()
        gaia_output_tmp['radius'] = radius_arr
        gaia_output_tmp['flux'] = flux_arr
        gaia_output_tmp['ls_'+band+'_bin'] = ls_mag_bins[mag_index]
        gaia_output.append(gaia_output_tmp)
        
        x, y, y_scatter = binned_stats(radius_arr.flatten(), flux_arr.flatten(), bins=np.linspace(0., np.floor(max_radius), nbins))
        radius_in_bin.append(x)
        flux_in_bin.append(y)
        flux_scatter_in_bin.append(y_scatter)
        
    gaia_output = vstack(gaia_output)

    ###################################################################################################

    for index in range(len(ls_mag_bins)):
        mask = gaia_output['ls_'+band+'_bin']==ls_mag_bins[index]
        plt.figure(figsize=(8, 6))
        plt.loglog(gaia_output['radius'][mask].T, gaia_output['flux'][mask].T, lw=1, alpha=0.1, c='C0')
        plt.plot(radius_in_bin[index], flux_in_bin[index], lw=1, alpha=1, c='C1', zorder=5)
        plt.axis([3, 300, .01, 500])
        plt.xlabel('Radius (arcsec)')
        plt.ylabel('SB (nmgy/sq.arcsec.)')
        plt.grid(alpha=0.5)
        plt.title('{} {} {}mag = {:.2f}'.format(field, output_name, band, ls_mag_bins[index]))
        plt.savefig('plots/more/{}_{}_{}mag_{:.2f}.png'.format(field, output_name, band, ls_mag_bins[index]))
        plt.close()

    plt.figure(figsize=(11, 8))
    for index in range(len(ls_mag_bins)):
        # normalize the flux to 13th magnitude stars
        norm = 10**((ls_mag_bins[index]-13)/2.5)
        plt.loglog(radius_in_bin[index], flux_in_bin[index]*norm, lw=1.5, alpha=1., 
                   label='{}mag = {:.2f}'.format(band, ls_mag_bins[index]), c='C'+str(index))
        # plt.errorbar(radius_in_bin[index], flux_in_bin[index]*norm, yerr=flux_scatter_in_bin[index],
        #              lw=1, alpha=1, c='C'+str(index), zorder=5)
    plt.title('{} {} {}-band'.format(field, output_name, band))
    plt.axis([.5, 100, .0002, 1000])
    plt.grid(alpha=0.5)
    plt.xlabel('Radius (arcsec)')
    plt.ylabel('SB (a.u.)')
    plt.legend()
    plt.savefig('plots/{}_{}_{}mag_average-bright.png'.format(field, output_name, band))
    plt.close()

    ###################################################################################################

    gaia_output.write('data/individual_profiles_{}_{}_{}-bright.fits'.format(field, output_name, band))

    t = Table()
    for index in range(len(flux_in_bin)):
        t['radius_{}_{:.2f}'.format(band, ls_mag_bins[index])] = radius_in_bin[index]
        t['flux_{}_{:.2f}'.format(band, ls_mag_bins[index])] = flux_in_bin[index]
        t['flux_scatter_{}_{:.2f}'.format(band, ls_mag_bins[index])] = flux_scatter_in_bin[index]
    t.write('data/average_profiles_{}_{}_{}-bright.txt'.format(field, output_name, band), format='ascii.commented_header')

    print()