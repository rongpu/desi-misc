from __future__ import division, print_function
import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from multiprocessing import Pool

tractor_fns = ['/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-xmm-combined.fits', '/global/cfs/cdirs/cosmo/work/legacysurvey/ibis/reductions/ibis-dr1/catalogs/tractor-cosmos-combined.fits']

tractor_cols = ['ibis_id_dr1', 'release', 'brickid', 'brickname', 'objid', 'maskbits', 'fitbits', 'type', 'ra', 'dec', 'ra_ivar', 'dec_ivar', 'dchisq', 'ebv', 'nearest_neighbor', 'ref_cat', 'ref_id', 'pmra', 'pmdec', 'parallax', 'pmra_ivar', 'pmdec_ivar', 'parallax_ivar', 'ref_epoch', 'gaia_phot_g_mean_mag', 'gaia_phot_g_mean_flux_over_error', 'gaia_phot_bp_mean_mag', 'gaia_phot_bp_mean_flux_over_error', 'gaia_phot_rp_mean_mag', 'gaia_phot_rp_mean_flux_over_error', 'gaia_astrometric_excess_noise', 'gaia_duplicated_source', 'gaia_phot_bp_rp_excess_factor', 'gaia_astrometric_sigma5d_max', 'gaia_astrometric_params_solved', 'flux_m411', 'flux_m438', 'flux_m464', 'flux_m490', 'flux_m517', 'flux_w1', 'flux_w2', 'flux_w3', 'flux_w4', 'flux_ivar_m411', 'flux_ivar_m438', 'flux_ivar_m464', 'flux_ivar_m490', 'flux_ivar_m517', 'flux_ivar_w1', 'flux_ivar_w2', 'flux_ivar_w3', 'flux_ivar_w4', 'fiberflux_m411', 'fiberflux_m438', 'fiberflux_m464', 'fiberflux_m490', 'fiberflux_m517', 'fibertotflux_m411', 'fibertotflux_m438', 'fibertotflux_m464', 'fibertotflux_m490', 'fibertotflux_m517', 'mw_transmission_w1', 'mw_transmission_w2', 'mw_transmission_w3', 'mw_transmission_w4', 'nobs_m411', 'nobs_m438', 'nobs_m464', 'nobs_m490', 'nobs_m517', 'nobs_w1', 'nobs_w2', 'nobs_w3', 'nobs_w4', 'rchisq_m411', 'rchisq_m438', 'rchisq_m464', 'rchisq_m490', 'rchisq_m517', 'rchisq_w1', 'rchisq_w2', 'rchisq_w3', 'rchisq_w4', 'fracflux_m411', 'fracflux_m438', 'fracflux_m464', 'fracflux_m490', 'fracflux_m517', 'fracflux_w1', 'fracflux_w2', 'fracflux_w3', 'fracflux_w4', 'fracmasked_m411', 'fracmasked_m438', 'fracmasked_m464', 'fracmasked_m490', 'fracmasked_m517', 'fracin_m411', 'fracin_m438', 'fracin_m464', 'fracin_m490', 'fracin_m517', 'anymask_m411', 'anymask_m438', 'anymask_m464', 'anymask_m490', 'anymask_m517', 'allmask_m411', 'allmask_m438', 'allmask_m464', 'allmask_m490', 'allmask_m517', 'wisemask_w1', 'wisemask_w2', 'psfsize_m411', 'psfsize_m438', 'psfsize_m464', 'psfsize_m490', 'psfsize_m517', 'psfdepth_m411', 'psfdepth_m438', 'psfdepth_m464', 'psfdepth_m490', 'psfdepth_m517', 'galdepth_m411', 'galdepth_m438', 'galdepth_m464', 'galdepth_m490', 'galdepth_m517', 'psfdepth_w1', 'psfdepth_w2', 'wise_coadd_id', 'sersic', 'sersic_ivar', 'shape_r', 'shape_r_ivar', 'shape_e1', 'shape_e1_ivar', 'shape_e2', 'shape_e2_ivar', 'lc_flux_w1', 'lc_flux_w2', 'lc_flux_ivar_w1', 'lc_flux_ivar_w2', 'lc_nobs_w1', 'lc_nobs_w2', 'lc_fracflux_w1', 'lc_fracflux_w2', 'lc_rchisq_w1', 'lc_rchisq_w2', 'lc_mjd_w1', 'lc_mjd_w2', 'lc_epoch_index_w1', 'lc_epoch_index_w2', 'brick_primary', 'bx', 'by', 'mjd_min', 'mjd_max', 'gaia_phot_g_n_obs', 'gaia_phot_bp_n_obs', 'gaia_phot_rp_n_obs', 'gaia_phot_variable_flag', 'gaia_astrometric_excess_noise_sig', 'gaia_astrometric_n_obs_al', 'gaia_astrometric_n_good_obs_al', 'gaia_astrometric_weight_al', 'gaia_a_g_val', 'gaia_e_bp_min_rp_val', 'apflux_m411', 'apflux_m438', 'apflux_m464', 'apflux_m490', 'apflux_m517', 'apflux_resid_m411', 'apflux_resid_m438', 'apflux_resid_m464', 'apflux_resid_m490', 'apflux_resid_m517', 'apflux_blobresid_m411', 'apflux_blobresid_m438', 'apflux_blobresid_m464', 'apflux_blobresid_m490', 'apflux_blobresid_m517', 'apflux_ivar_m411', 'apflux_ivar_m438', 'apflux_ivar_m464', 'apflux_ivar_m490', 'apflux_ivar_m517', 'apflux_masked_m411', 'apflux_masked_m438', 'apflux_masked_m464', 'apflux_masked_m490', 'apflux_masked_m517', 'apflux_w1', 'apflux_w2', 'apflux_w3', 'apflux_w4', 'apflux_resid_w1', 'apflux_resid_w2', 'apflux_resid_w3', 'apflux_resid_w4', 'apflux_ivar_w1', 'apflux_ivar_w2', 'apflux_ivar_w3', 'apflux_ivar_w4', 'ngood_m411', 'ngood_m438', 'ngood_m464', 'ngood_m490', 'ngood_m517', 'nea_m411', 'nea_m438', 'nea_m464', 'nea_m490', 'nea_m517', 'blob_nea_m411', 'blob_nea_m438', 'blob_nea_m464', 'blob_nea_m490', 'blob_nea_m517', 'psfdepth_w3', 'psfdepth_w4', 'wise_x', 'wise_y']

def read_hsc_catalog(brickname):
    hsc_fn = '/pscratch/sd/r/rongpu/ibis-tractor/ibis-dr1-hsc-wide-forced/forced-brick/{}/tractor-forced-{}.fits'.format(brickname[:3], brickname)
    if os.path.isfile(hsc_fn):
        hsc = Table(fitsio.read(hsc_fn))
    else:
        print(hsc_fn, 'does not exist!')
        return None

    dup_columns = np.array(hsc.colnames)[np.in1d(hsc.colnames, tractor_cols)]
    hsc.remove_columns(dup_columns)
    hsc.rename_column('ls_id_dr11', 'ibis_id_dr1')

    return hsc

for tractor_fn in tractor_fns:
    cat = Table(fitsio.read(tractor_fn, columns=['ibis_id_dr1', 'brickname']))
    print(len(cat))

    bricknames = np.unique(cat['brickname'])    
    cat.remove_column('brickname')

    n_processes = 128
    with Pool(processes=n_processes) as pool:
        res = pool.map(read_hsc_catalog, bricknames)

    # Remove None elements from the list
    for index in range(len(res)-1, -1, -1):
        if res[index] is None:
            res.pop(index)

    hsc = vstack(res)

    ibis_id_dr1_original = cat['ibis_id_dr1'].copy()
    cat = join(cat, hsc, keys='ibis_id_dr1', join_type='left').filled(0)
    if not np.all(cat['ibis_id_dr1']==ibis_id_dr1_original):
        if len(cat)!=len(ibis_id_dr1_original) or not np.all(np.unique(cat['ibis_id_dr1'])==np.unique(ibis_id_dr1_original)):
            raise ValueError
        reverse_sort = np.array(ibis_id_dr1_original).argsort().argsort()
        cat = cat[np.argsort(cat['ibis_id_dr1'])[reverse_sort]]
    assert np.all(cat['ibis_id_dr1']==ibis_id_dr1_original)

    cat.write(tractor_fn.replace('-combined.fits', '-hsc-wide-forced.fits'))
