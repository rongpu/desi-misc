from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool


sweep_columns = ['ibis_id_dr1', 'release', 'brickid', 'brickname', 'objid', 'maskbits', 'fitbits', 'type', 'ra', 'dec', 'ra_ivar', 'dec_ivar', 'dchisq', 'ebv', 'nearest_neighbor', 'ref_cat', 'ref_id', 'pmra', 'pmdec', 'parallax', 'pmra_ivar', 'pmdec_ivar', 'parallax_ivar', 'ref_epoch', 'gaia_phot_g_mean_mag', 'gaia_phot_g_mean_flux_over_error', 'gaia_phot_bp_mean_mag', 'gaia_phot_bp_mean_flux_over_error', 'gaia_phot_rp_mean_mag', 'gaia_phot_rp_mean_flux_over_error', 'gaia_astrometric_excess_noise', 'gaia_duplicated_source', 'gaia_phot_bp_rp_excess_factor', 'gaia_astrometric_sigma5d_max', 'gaia_astrometric_params_solved', 'flux_M411', 'flux_M438', 'flux_M464', 'flux_M490', 'flux_M517', 'flux_w1', 'flux_w2', 'flux_w3', 'flux_w4', 'flux_ivar_M411', 'flux_ivar_M438', 'flux_ivar_M464', 'flux_ivar_M490', 'flux_ivar_M517', 'flux_ivar_w1', 'flux_ivar_w2', 'flux_ivar_w3', 'flux_ivar_w4', 'fiberflux_M411', 'fiberflux_M438', 'fiberflux_M464', 'fiberflux_M490', 'fiberflux_M517', 'fibertotflux_M411', 'fibertotflux_M438', 'fibertotflux_M464', 'fibertotflux_M490', 'fibertotflux_M517', 'mw_transmission_w1', 'mw_transmission_w2', 'mw_transmission_w3', 'mw_transmission_w4', 'nobs_M411', 'nobs_M438', 'nobs_M464', 'nobs_M490', 'nobs_M517', 'nobs_w1', 'nobs_w2', 'nobs_w3', 'nobs_w4', 'rchisq_M411', 'rchisq_M438', 'rchisq_M464', 'rchisq_M490', 'rchisq_M517', 'rchisq_w1', 'rchisq_w2', 'rchisq_w3', 'rchisq_w4', 'fracflux_M411', 'fracflux_M438', 'fracflux_M464', 'fracflux_M490', 'fracflux_M517', 'fracflux_w1', 'fracflux_w2', 'fracflux_w3', 'fracflux_w4', 'fracmasked_M411', 'fracmasked_M438', 'fracmasked_M464', 'fracmasked_M490', 'fracmasked_M517', 'fracin_M411', 'fracin_M438', 'fracin_M464', 'fracin_M490', 'fracin_M517', 'anymask_M411', 'anymask_M438', 'anymask_M464', 'anymask_M490', 'anymask_M517', 'allmask_M411', 'allmask_M438', 'allmask_M464', 'allmask_M490', 'allmask_M517', 'wisemask_w1', 'wisemask_w2', 'psfsize_M411', 'psfsize_M438', 'psfsize_M464', 'psfsize_M490', 'psfsize_M517', 'psfdepth_M411', 'psfdepth_M438', 'psfdepth_M464', 'psfdepth_M490', 'psfdepth_M517', 'galdepth_M411', 'galdepth_M438', 'galdepth_M464', 'galdepth_M490', 'galdepth_M517', 'psfdepth_w1', 'psfdepth_w2', 'wise_coadd_id', 'sersic', 'sersic_ivar', 'shape_r', 'shape_r_ivar', 'shape_e1', 'shape_e1_ivar', 'shape_e2', 'shape_e2_ivar']
lc_columns = ['ibis_id_dr1', 'release', 'brickid', 'objid', 'lc_flux_w1', 'lc_flux_w2', 'lc_flux_ivar_w1', 'lc_flux_ivar_w2', 'lc_nobs_w1', 'lc_nobs_w2', 'lc_fracflux_w1', 'lc_fracflux_w2', 'lc_rchisq_w1', 'lc_rchisq_w2', 'lc_mjd_w1', 'lc_mjd_w2', 'lc_epoch_index_w1', 'lc_epoch_index_w2']
extra_columns = ['ibis_id_dr1', 'release', 'brickid', 'objid', 'brick_primary', 'bx', 'by', 'mjd_min', 'mjd_max', 'gaia_phot_g_n_obs', 'gaia_phot_bp_n_obs', 'gaia_phot_rp_n_obs', 'gaia_phot_variable_flag', 'gaia_astrometric_excess_noise_sig', 'gaia_astrometric_n_obs_al', 'gaia_astrometric_n_good_obs_al', 'gaia_astrometric_weight_al', 'gaia_a_g_val', 'gaia_e_bp_min_rp_val', 'apflux_M411', 'apflux_M438', 'apflux_M464', 'apflux_M490', 'apflux_M517', 'apflux_resid_M411', 'apflux_resid_M438', 'apflux_resid_M464', 'apflux_resid_M490', 'apflux_resid_M517', 'apflux_blobresid_M411', 'apflux_blobresid_M438', 'apflux_blobresid_M464', 'apflux_blobresid_M490', 'apflux_blobresid_M517', 'apflux_ivar_M411', 'apflux_ivar_M438', 'apflux_ivar_M464', 'apflux_ivar_M490', 'apflux_ivar_M517', 'apflux_masked_M411', 'apflux_masked_M438', 'apflux_masked_M464', 'apflux_masked_M490', 'apflux_masked_M517', 'apflux_w1', 'apflux_w2', 'apflux_w3', 'apflux_w4', 'apflux_resid_w1', 'apflux_resid_w2', 'apflux_resid_w3', 'apflux_resid_w4', 'apflux_ivar_w1', 'apflux_ivar_w2', 'apflux_ivar_w3', 'apflux_ivar_w4', 'ngood_M411', 'ngood_M438', 'ngood_M464', 'ngood_M490', 'ngood_M517', 'nea_M411', 'nea_M438', 'nea_M464', 'nea_M490', 'nea_M517', 'blob_nea_M411', 'blob_nea_M438', 'blob_nea_M464', 'blob_nea_M490', 'blob_nea_M517', 'psfdepth_w3', 'psfdepth_w4', 'wise_x', 'wise_y']

fns = sorted(glob.glob('/dvs_ro/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/tractor/*/tractor-*.fits'))
print(len(fns))

def read_catalog(index):
    fn = fns[index]
    cat = Table(fitsio.read(fn))
    cat = cat[cat['brick_primary']]
    return cat

n_processes = 128
with Pool(processes=n_processes) as pool:
    res = pool.map(read_catalog, np.arange(len(fns)))
cat = vstack(res)
print(len(cat))

cat.rename_column('ls_id_dr11', 'ibis_id_dr1')

cat_all = cat.copy()

mask = cat_all['ra']<100
cat = cat_all[mask].copy()
print('XMM', len(cat))
tractor = cat[sweep_columns].copy()
tractor_lc = cat[lc_columns].copy()
tractor_extra = cat[extra_columns].copy()
tractor.write('/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-xmm-combined.fits')
tractor_lc.write('/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-xmm-combined-lc.fits')
tractor_extra.write('/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-xmm-combined-extra.fits')

mask = cat_all['ra']>100
cat = cat_all[mask].copy()
print('COSMOS', len(cat))
tractor = cat[sweep_columns].copy()
tractor_lc = cat[lc_columns].copy()
tractor_extra = cat[extra_columns].copy()
tractor.write('/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-cosmos-combined.fits')
tractor_lc.write('/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-cosmos-combined-lc.fits')
tractor_extra.write('/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-cosmos-combined-extra.fits')
