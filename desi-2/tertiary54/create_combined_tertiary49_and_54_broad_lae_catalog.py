# Create the combined broadly-selected LAE catalog from tertiary49 and teritary54 for studying the LAE selection efficiencies.

import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp


params = {'legend.fontsize': 'medium',
          'axes.labelsize': 'medium',
          'axes.titlesize': 'medium',
          'xtick.labelsize': 'medium',
          'ytick.labelsize': 'medium',
          'figure.facecolor': 'w'}
plt.rcParams.update(params)


sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import match_coord


max_radius = 1.5  # deg

lyman_alpha_z = {'M411': [2.3, 2.51], 'M438': [2.51, 2.72], 'M464': [2.72, 2.93], 'M490': [2.93, 3.15], 'M517': [3.15, 3.36]}
ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']

ext_coeffs = {
'g': 3.240, 'r': 2.276, 'i': 1.633, 'z': 1.263, 'y': 1.076,
'M411': 4.290, 'M438': 4.099, 'M464': 3.877, 'M490': 3.634, 'M517': 3.389}

# compute magnitudes
fnodet = 0.001
magnodet = 22.5 - 2.5 * np.log10(fnodet) # 30

redshift_columns = ['Z','ZERR','ZWARN','CHI2','SPECTYPE','DELTACHI2','PETAL_LOC','DEVICE_LOC','LOCATION','FIBER','COADD_FIBERSTATUS','TARGET_RA','TARGET_DEC','OBJTYPE','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY','TILEID','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT','COADD_NUMTILE','TSNR2_BGS_B','TSNR2_ELG_B','TSNR2_LRG_B','TSNR2_QSO_B','TSNR2_BGS','TSNR2_ELG','TSNR2_LRG','TSNR2_QSO','rr_fn','lya_flux','lya_flux_err','coadd_fn','EFFTIME_LRG','low_z','lae','Z_daily','ZWARN_daily','DELTACHI2_daily','SPECTYPE_daily','ra','dec','CENTER_OFFSET','M411_sel','M438_sel','M464_sel','M490_sel','M517_sel','lae_sel']


# ---------
# # tertiary49

targets = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/desi2/xmm_highz/tertiary49_xmm_lae_targets.fits', columns=['FIELD_RADIUS_OFFSET']))
mask = targets['FIELD_RADIUS_OFFSET']<max_radius
target_density_49 = np.sum(mask) / (np.pi * max_radius**2)
print('Tertiery49 target density: {:.1f} per sq deg\n'.format(target_density_49))

cat = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/catalogs/tertiary49_ibis_lae_truth-20260127-add_target_info-add_lyaflux.fits'))
print(len(cat), len(np.unique(cat['TARGETID'])))

cat.rename_columns(['FIELD_RADIUS_OFFSET', 'RA', 'DEC'], ['CENTER_OFFSET', 'ra', 'dec'])
mask = cat['CENTER_OFFSET']<max_radius
cat = cat[mask]
print(len(cat))

tertiary49 = cat[redshift_columns].copy()


# ---------
# # tertiary54

primary_targ = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/tertiary54_lae_targets.fits'))
print(len(primary_targ))
mask = primary_targ['CENTER_OFFSET']<max_radius
primary_targ = primary_targ[mask]
print(len(primary_targ))

pp = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/catalogs/tertiary54_ibis_primary_lae_truth-20260413.fits'))
print(len(pp), len(np.unique(pp['TARGETID'])))
mask = pp['CENTER_OFFSET']<max_radius
pp = pp[mask]
print(len(pp))

filler_targ = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/tertiary54_lae_filler_targets.fits'))
print(len(filler_targ))
mask = filler_targ['CENTER_OFFSET']<max_radius
filler_targ = filler_targ[mask]
print(len(filler_targ))

ff = Table(fitsio.read('/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/catalogs/tertiary54_ibis_filler_lae_truth-20260413.fits'))
print(len(ff), len(np.unique(ff['TARGETID'])))
mask = ff['CENTER_OFFSET']<max_radius
ff = ff[mask]
print(len(ff))

target_density_54 = (len(primary_targ) + len(filler_targ)) / (np.pi * max_radius**2)
print('Tertiery54 target density: {:.1f} per sq deg\n'.format(target_density_54))

# Rebalance the sample by downsampling the primary sample
primary_completeness = len(pp)/len(primary_targ)
filler_completeness = len(ff)/len(filler_targ)
print('Primary completeness: {:.2f}%'.format(100*primary_completeness))
print('Filler completeness: {:.2f}%'.format(100*filler_completeness))
print()

n = int(np.round(len(pp) / (primary_completeness / filler_completeness)))
np.random.seed(666)
idx = np.sort(np.random.choice(len(pp), size=n, replace=False))
pp = pp[idx]

primary_completeness = len(pp)/len(primary_targ)
filler_completeness = len(ff)/len(filler_targ)
print('Primary completeness: {:.2f}%'.format(100*primary_completeness))
print('Filler completeness: {:.2f}%'.format(100*filler_completeness))

tertiary54 = vstack([pp, ff], join_type='inner')
print(len(tertiary54))

tertiary54 = tertiary54[redshift_columns]


# # Combined

tertiary49['tertiary'] = 49
tertiary54['tertiary'] = 54

cat = vstack([tertiary49, tertiary54], join_type='exact')
print(len(cat), np.sum(cat['lae']))

effective_area = np.sum(cat['tertiary']==49)/target_density_49 \
               + np.sum(cat['tertiary']==54)/target_density_54
print('Effective area: {:.4f} sq deg\n'.format(effective_area))

cat_all = cat.copy()


photom_dirs = [
    '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_fields/',
    '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/subset-25.0/',
    '/dvs_ro/cfs/cdirs/desicollab/users/rongpu/data/ibis/deep_field_subsets/subset-24.8/',
]

for photom_dir in photom_dirs:
    
    cat = cat_all.copy()

    print(photom_dir)
    ibis_stack = []
    for field in ['cosmos', 'xmm']:
        ibis_fn = os.path.join(photom_dir, 'tractor_{}.fits'.format(field))
        hsc_fn = os.path.join(photom_dir, 'tractor_{}-hsc_wide.fits'.format(field))
        ibis = Table(fitsio.read(ibis_fn))
        hsc = Table(fitsio.read(hsc_fn))
        
        assert np.all(ibis['ibis_id']==hsc['ibis_id'])
        # print(len(ibis), len(hsc))
        # print(np.intersect1d(ibis.colnames, hsc.colnames))
        hsc.remove_columns(['ibis_id', 'ra', 'dec'])
        ibis = hstack([ibis, hsc])
        # print(len(ibis))

        ibis_stack.append(ibis)

    ibis = vstack(ibis_stack)
    
    print('Cross-matching')
    idx1, idx2, d2d, d_ra, d_dec = match_coord(ibis['ra'], ibis['dec'], cat['ra'], cat['dec'], search_radius=0.5, plot_q=True)
    cross_matching_completeness = len(idx2)/len(cat)
    print(len(idx2))
    print('Cross-matching completeness: {:.2f}%'.format(100*cross_matching_completeness))
    ibis = ibis[idx1]
    cat = cat[idx2]
    print(len(cat))
    
    print(np.intersect1d(ibis.colnames, cat.colnames))
    for col in np.intersect1d(ibis.colnames, cat.colnames):
        ibis.remove_column(col)
    
    cat = hstack([cat, ibis])
    
    # clip the magnitudes at 30.0
    ibis_filters = ['M411', 'M438', 'M464', 'M490', 'M517']
    for band in ibis_filters:
        cat[band+'mag'] = 22.5 - 2.5*np.log10(np.clip(cat['flux_{}'.format(band)]*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))
        cat[band+'fibmag'] = 22.5 - 2.5*np.log10(np.clip(cat['fiberflux_{}'.format(band)]*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))
    for band in ['g', 'r', 'i', 'z', 'y']:
        mask = ~np.isnan(cat['{}_cmodel_flux'.format(band)])
        cat[band+'mag'] = -99.
        cat[band+'mag'][mask] = 22.5 - 2.5*np.log10(np.clip(cat['{}_cmodel_flux'.format(band)]/3630.78*10**(0.4*ext_coeffs[band]*cat['ebv']), fnodet, None))[mask]
        cat[band+'mag'][~mask] = 30.0

    # leave-one-out average magnitude
    for band in ibis_filters:
        cat[band+'_ref_flux'] = 0.
        for band1 in ibis_filters:
            if band1!=band:
                cat[band+'_ref_flux'] += cat['flux_'+band1] * 10**(0.4*ext_coeffs[band1]*cat['ebv']) / 4
        cat[band+'_ref_mag'] = 22.5 - 2.5*np.log10(np.clip(cat[band+'_ref_flux'], fnodet, None))
    
    if photom_dir.split('/')[-2]=='deep_fields':
        output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/tertiary49_and_54/tertiary49_and_54_broad_lae_truth-deep.fits'
    else:
        output_fn = '/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary54/tertiary49_and_54/tertiary49_and_54_broad_lae_truth-{}.fits'.format(photom_dir.split('/')[-2])
    cat.write(output_fn, overwrite=False)
