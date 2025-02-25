from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
from astropy.io import fits

from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
from matplotlib.colors import LogNorm
import healpy as hp

from multiprocessing import Pool


nmad = lambda x: 1.4826 * np.median(np.abs(x-np.median(x)))

pixscale = 0.262  # arcsec

psfex_dir = '/global/cfs/cdirs/cosmo/work/legacysurvey/dr11/calib/psfex'
gaia_dir = '/global/cfs/cdirs/cosmo/data/gaia/dr3/healpix'

ccd = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/dr11/survey-ccds-decam-dr11-psfex-stats.fits'))
print(len(ccd))

mask = ccd['ccd_cuts']==0
ccd = ccd[mask]
print(len(ccd))

mask = ccd['dec']>-28
ccd = ccd[mask]
print(len(ccd))

#################################
mask = ccd['fwhm_psfex']<3.
ccd = ccd[mask]
print(len(ccd))
#################################

ccd['r_cent_recentered'] = np.sqrt((ccd['x_cent']-np.median(ccd['x_cent']))**2 + (ccd['y_cent']-np.median(ccd['y_cent']))**2)

# with open('/global/u2/r/rongpu/moregit/legacypipe/py/legacyzpts/data/decam-bad_expid.txt') as f:
#     texts = list(map(str.strip, f.readlines()))
# known_bad_expnum_list = []
# for text in texts:
#     if len(text)>0 and (text[0]!='#'):
#         known_bad_expnum_list.append(int(text.replace('-', ' ').split()[0]))
# print('bad exposure list:', len(known_bad_expnum_list))

# mask = np.in1d(ccd['expnum'], known_bad_expnum_list)
# print('in bad exposure list:', np.sum(mask))
# ccd = ccd[~mask]


def plot_psf(ccd_index, show=False):

    expnum = ccd['expnum'][ccd_index]

    # plot_path = '/global/cfs/cdirs/cosmo/www/temp/rongpu/dr11/weird_psfs/decam_{}.png'.format(expnum)
    plot_path = '/global/cfs/cdirs/cosmo/www/temp/rongpu/dr11/bad_psfs/decam_{}.png'.format(expnum)
    if os.path.isfile(plot_path) and show==False:
        return None

    ############################ Load image ############################

    image_dir = '/global/cfs/cdirs/cosmo/staging/'
    fn = os.path.join(image_dir, ccd['image_filename'][ccd_index].strip())
    hdulist = fits.open(fn)
    w = wcs.WCS(hdulist[ccd['image_hdu'][ccd_index]].header)
    img = hdulist[ccd['image_hdu'][ccd_index]].data

    # naive sky subtraction
    mask = (img<np.percentile(img.flatten(), 85))
    img = img - np.median(img[mask].flatten())

    sky_nmad = nmad(img.flatten())

    ############################ Find Gaia stars ############################
    
    gaia_nside = 32
    gaia_npix = hp.nside2npix(gaia_nside)
    gaia_hp_ra, gaia_hp_dec = hp.pix2ang(gaia_nside, np.arange(gaia_npix), nest=True, lonlat=True)

    sky1 = SkyCoord(ccd['ra'][[ccd_index]]*u.degree, ccd['dec'][[ccd_index]]*u.degree, frame='icrs')
    sky2 = SkyCoord(gaia_hp_ra*u.degree, gaia_hp_dec*u.degree, frame='icrs')
    sep = sky1.separation(sky2)
    sep = sep.to_value('deg')
    gaia_hp_index = np.argmin(sep)

    gaia_fn = (5-len(str(gaia_hp_index)))*'0'+str(gaia_hp_index)
    gaia = Table.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)))

    # Gaia magnitude cuts
    mask = (gaia['PHOT_G_MEAN_MAG']>13.) & (gaia['PHOT_G_MEAN_MAG']<19)
    gaia = gaia[mask]

    # Select stars within the CCD
    gaia['x'], gaia['y'] = w.wcs_world2pix(gaia['RA'], gaia['DEC'], True)
    mask = ~((gaia['x']<31) | (gaia['y']<31) | (gaia['x']>img.shape[1]-31) | (gaia['y']>img.shape[0]-31))
    gaia = gaia[mask]
    # print('{} gaia stars'.format(len(gaia)))

    # shuffle
    np.random.seed(expnum)
    idx = np.random.choice(len(gaia), size=len(gaia), replace=False)
    gaia = gaia[idx]

    # Select the brightest Gaia stars
    if len(gaia)>7:
        idx = np.argsort(gaia['PHOT_G_MEAN_MAG'])
        idx1 = idx[:3]
        np.random.seed(ccd_index)
        idx2 = np.random.choice(idx[3:], size=4, replace=False)
        idx = np.concatenate([idx1, idx2])
        gaia = gaia[idx]

    ############################ Load PSFEx ############################

    image_filename = ccd['image_filename'][ccd_index]
    psfex_filename = image_filename.replace('.fits.fz', '-psfex.fits')
    psfex_path = os.path.join(psfex_dir, psfex_filename)
    data = Table(fitsio.read(psfex_path, ext=1))
    psf = np.array(data['psf_mask'][0][0])

    ############################ Make plot ############################

    y_size, x_size = psf.shape
    fig, ax = plt.subplots(3, 3, figsize=(12, 12))
    ax_flat = ax.flatten()
    if psf.max()>0:
        ax_flat[0].imshow(psf, cmap='viridis', norm=LogNorm(vmin=1e-3*psf.max(), vmax=psf.max()), origin='lower')
        ax_flat[0].axvline((x_size-1)/2, color='0.7', lw=1, ls='--')
        ax_flat[0].axhline((y_size-1)/2, color='0.7', lw=1, ls='--')
        ax_flat[0].text(0, 0, 'log scale', color='0.7', fontsize=13)
        ax_flat[1].imshow(psf, cmap='gray_r', origin='lower')
        ax_flat[1].axvline((x_size-1)/2, color='0.7', lw=1, ls='--')
        ax_flat[1].axhline((y_size-1)/2, color='0.7', lw=1, ls='--')
    for index in range(0, len(gaia)):
        x, y = int(gaia['x'][index]), int(gaia['y'][index])
        tmp = img[y-30:y+30, x-30:x+30]
        # ax_flat[index+2].imshow(tmp, norm=LogNorm(vmin=5*sky_nmad), origin='lower')
        ax_flat[index+2].imshow(tmp, origin='lower', cmap='gray_r')
        ax_flat[index+2].axis('off')
    title_text = ''
    # if expnum in known_bad_expnum_list:
    #     title_text += 'known_bad '
    title_text += '{} band; ccd_cuts={}'.format(ccd['filter'][ccd_index], ccd['ccd_cuts'][ccd_index])
    title_text += '\n fwhm={:.2f} rmax={:.2f} rcenter={:.1f}'.format(ccd['fwhm_psfex'][ccd_index], ccd['r_max'][ccd_index], ccd['r_cent_recentered'][ccd_index])
    ax_flat[0].set_title(title_text)
    ax_flat[1].set_title('expnum={}\npropid={}'.format(ccd['expnum'][ccd_index], ccd['propid'][ccd_index]))
    text = 'ra={:.2f} dec={:.2f}'.format(ccd['ra'][ccd_index], ccd['dec'][ccd_index])
    text += '\n ellipticity={:.2f} centroid_max_ratio={:.3f}'.format(ccd['ellipticity'][ccd_index], ccd['centroid_max_ratio'][ccd_index])
    ax_flat[2].set_title(text)
    plt.tight_layout()
    if show==False:
        plt.savefig(plot_path)
        plt.close()
    else:
        plt.show()

    return None


# plot_psf(0, show=True)

# mask = (ccd['r_cent_recentered']>0.6)
# print(np.sum(mask), '{:.3f}%'.format(100*np.sum(mask)/len(mask)))
# idx = np.where(mask)[0]


with open('/global/u2/r/rongpu/temp/decam-bad_expid-more.txt') as f:
    texts = list(map(str.strip, f.readlines()))
known_bad_expnum_list = []
for text in texts:
    if len(text)>0 and (text[0]!='#'):
        known_bad_expnum_list.append(int(text.replace('-', ' ').split()[0]))
print('bad exposure list:', len(known_bad_expnum_list))
mask = np.in1d(ccd['expnum'], known_bad_expnum_list)
print('in bad exposure list:', np.sum(mask))
ccd = ccd[mask]

n_process = 128
with Pool(processes=n_process) as pool:
    res = pool.map(plot_psf, np.arange(len(ccd)))
