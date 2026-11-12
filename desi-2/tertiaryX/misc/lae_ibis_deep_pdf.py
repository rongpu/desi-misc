# Created by Claude
# LAE selection diagnostic plots: one PDF page per band, two rows per page
# (top: tertiary49 truth-table selection; bottom: tertiaryX catalog selection)

import sys, os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from astropy.table import Table
import fitsio
from scipy.stats import norm
from astropy.coordinates import SkyCoord
from astropy import units as u

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord, match_self

params = {'legend.fontsize': 'medium',
          'axes.labelsize': 'medium',
          'axes.titlesize': 'medium',
          'xtick.labelsize': 'medium',
          'ytick.labelsize': 'medium',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)


def flux_threshold(f_b, depth, c0, prob_tail, sigma_c=0., rescale_factor=1, include_color_lim=False):
    '''
    Continuum flux threshold f_ref_lim for selecting narrow-band excess objects
    (copied unchanged from the original script; see there for the full docstring).
    '''
    f_b = np.asarray(f_b, dtype=float)
    z = norm.isf(prob_tail)  # gaussian sigma at tail probability prob_tail
    a = 10**(0.4 * c0)
    sigma_b = rescale_factor * 10**(-0.4 * (depth - 22.5)) / 5  # error on selection-band flux
    sigma_x = sigma_b * np.sqrt(1 + a**2 / 4)  # error on f_b - a*f_ref
    # solve f_b - y = z * sqrt(sigma_x**2 + (k*y)**2) for y = a * f_ref_lim (root with y < f_b)
    k = (10**(0.4 * z * sigma_c) - 1) / z
    A = 1 - (z * k)**2
    if A <= 0:
        raise ValueError('sigma_c too large for this prob_tail')
    disc = f_b**2 - A * (f_b**2 - (z * sigma_x)**2)
    f_ref_lim = (f_b - np.sqrt(disc)) / A / a

    if not include_color_lim:
        return f_ref_lim
    else:
        with np.errstate(divide='ignore', invalid='ignore'):
            color_lim = np.where(f_ref_lim > 0, 2.5 * np.log10(f_b / f_ref_lim), np.inf)
        color_lim = np.where(f_b > 0, color_lim, np.nan)
        return f_ref_lim, color_lim


median_depths = {'M411': 25.39, 'M438': 25.45, 'M464': 25.53, 'M490': 25.54, 'M517': 25.59}
contam_median_colors = {'M411': -0.07, 'M438': -0.02, 'M464': 0.03, 'M490': 0.1, 'M517': 0.12}

cat = Table(fitsio.read('/pscratch/sd/r/rongpu/tmp/tertiaryX/tertiaryX_lae_targets-all.fits'))
print(len(cat))
mask = cat['CENTER_OFFSET'] < 1.45
cat = cat[mask]
print(len(cat))

tt = Table(fitsio.read('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/catalogs/tertiary49_ibis_lae_truth-20260127-add_target_info-add_lyaflux-newcol.fits'))
print(len(tt))

field_ra, field_dec, field_rad_deg = 35.75, -4.75, 1.8

# cut on radius
field_c = SkyCoord(ra=field_ra * u.deg, dec=field_dec * u.deg, frame='icrs')
cs = SkyCoord(ra=tt['RA'] * u.deg, dec=tt['DEC'] * u.deg, frame='icrs')
tt['CENTER_OFFSET'] = cs.separation(field_c).deg
mask = tt['CENTER_OFFSET'] < 1.45
tt = tt[mask]

################################### Remove duplicates ###################################
score = tt['EFFTIME_LRG'].copy()
n_duplicates = 1

while n_duplicates>0:
    n_duplicates, idx1, idx2 = match_self(tt['TARGET_RA'], tt['TARGET_DEC'], search_radius=1., return_indices=True, plot_q=False)
    if len(idx1)==0:
        continue
    mask_remove1 = (score[idx1] < score[idx2]) | (score[idx1]==0)
    mask_remove2 = (score[idx2] < score[idx1]) | (score[idx2]==0)
    idx_remove = np.unique(np.concatenate([idx1[mask_remove1], idx2[mask_remove2]]))
    print(len(idx_remove), len(idx_remove)/len(idx1))
    mask = np.full(len(tt), True)
    mask[idx_remove] = False
    print(np.sum(mask), np.sum(mask)/len(mask))
    tt = tt[mask]
    score = score[mask]

print()
print(len(tt))

idx1, idx2, d2d, d_ra, d_dec = match_coord(tt['RA'], tt['DEC'], cat['ra'], cat['dec'], search_radius=1., plot_q=False)
print(len(idx1)/len(tt))
tt = tt[idx1]
cat = cat[idx2]
print(len(cat))

good = tt['lae'].copy()

PROB_TAIL = 0.0001

LAE_CONFIGS = {
    'M411': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': -0.075, 'color_offset_sigma': 0.087},
    'M438': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': -0.085, 'color_offset_sigma': 0.08},
    'M464': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': -0.125, 'color_offset_sigma': 0.098},
    'M490': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': -0.03, 'color_offset_sigma': 0.121},
    'M517': {'gr_colmax': 0.7, 'ri_colmax': 0.35, 'color_offset': 0., 'color_offset_sigma': 0.178},
}

# (row title, table holding the {band}_sel column, style of the "not confirmed" points)
rows = [
    ('tertiary49 selection', tt, dict(marker='.', ms=3, alpha=0.6)),
    ('new tertiary selection', cat, dict(marker='o', ms=3, alpha=0.8)),
]

with PdfPages('/global/cfs/cdirs/desicollab/users/rongpu/tmp/claude_dump/LAE_IBIS_deep.pdf') as pdf:
    for band, cfg in LAE_CONFIGS.items():
        print(band)

        # Color-limit curve (same depth/contaminant color/scatter as in select_lae)
        mag_list = np.linspace(21, 26, 100)
        flux_list = 10**(-0.4 * (mag_list - 22.5))
        _, color_lim = flux_threshold(flux_list, depth=median_depths[band], c0=contam_median_colors[band],
                                      prob_tail=PROB_TAIL, sigma_c=cfg['color_offset_sigma'], include_color_lim=True)
        color_lim += cfg['color_offset']

        color = cat[f'{band}_ref_mag'] - cat[f'{band}mag']

        fig, axes = plt.subplots(2, 2, figsize=(12, 12))
        for ax, (title, sel_table, bg) in zip(axes, rows):
            mask1 = sel_table[band+'_sel'] & good
            mask2 = sel_table[band+'_sel'] & (~good)
            print(title, np.sum(mask1), np.sum(mask2), np.sum(mask1)/np.sum(mask1 | mask2))
            n1, n2 = np.sum(mask1), np.sum(mask2)
            label1 = f'confirmed LAE {100*n1/(n1+n2):.1f}% ({n1}/{n1+n2})'
            label2 = f'not confirmed LAE {100*n2/(n1+n2):.1f}% ({n2}/{n1+n2})'

            ax[0].plot(cat[f'{band}mag'][mask1], color[mask1], 'o', ms=2.7, mfc='C3', mec='w', mew=0.3, alpha=0.9, label=label1, zorder=2)
            ax[0].plot(cat[f'{band}mag'][mask2], color[mask2], ls='', mew=0, label=label2, color='0.6', zorder=0, **bg)
            ax[0].plot(mag_list, color_lim, color='k', lw=1.5, zorder=3)
            ax[0].set_xlabel(f'{band}mag')
            ax[0].set_ylabel(f'{band}_ref_mag - {band}mag')
            ax[0].axis([21.5, 25.25, -1, 3.5])
            ax[0].legend(markerscale=2, loc='upper left')
            ax[0].grid(alpha=0.3)
            ax[0].set_title(f'{band}: {title}')
            ax[1].plot((cat['rmag'] - cat['imag'])[mask1], (cat['gmag'] - cat['rmag'])[mask1], 'o', ms=2.7, mfc='C3', mec='w', mew=0.3, alpha=0.9, zorder=2)
            ax[1].plot((cat['rmag'] - cat['imag'])[mask2], (cat['gmag'] - cat['rmag'])[mask2], ls='', mew=0, color='0.6', zorder=0, **bg)
            ax[1].axvline(cfg['ri_colmax'], color='k', ls='--', lw=1.2, zorder=3)
            ax[1].axhline(cfg['gr_colmax'], color='k', ls='--', lw=1.2, zorder=3)
            ax[1].set_xlabel('r - i')
            ax[1].set_ylabel('g - r')
            ax[1].axis([-1., 1., -1, 1])
            ax[1].grid(alpha=0.3)
            ax[1].set_title(f'{band}: {title}')
        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
