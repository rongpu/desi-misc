import sys, os, glob, time, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits
import healpy as hp


fn = '/dvs_ro/cfs/cdirs/desi/spectro/redux/matterhorn/zcatalog/v2/zall/zall-pix-matterhorn.fits'
fn_extra = fn.replace('.fits', '-extra.fits')

cat = Table(fitsio.read(fn))
cat_extra = Table(fitsio.read(fn_extra))
assert np.all(cat['TARGETID']==cat_extra['TARGETID'])

for col in cat_extra.colnames:
    if col in cat.colnames:
        cat_extra.remove_column(col)
cat = hstack([cat, cat_extra])
print(len(cat))



import astropy
from desitarget.targetmask import desi_mask

sys.path.append(os.path.expanduser('~/moregit/desispec/py/desispec'))
import validredshifts


cat['id'] = np.arange(len(cat))

cat_stack = []

for survey in np.unique(cat['SURVEY']):

    print(survey)
    
    mask = cat['SURVEY']==survey
    zcat = cat[mask].copy()

    # LSS redshift quality cuts
    if survey=='main':
        zqual = validredshifts.actually_validate(zcat, populate_missing_columns=True)
    else:
        zqual = validredshifts.actually_validate(zcat, ignore_lya=True, populate_missing_columns=True)

    # GOOD_SPEC: true if it is a science spectrum with good hardware status
    good_spec = validredshifts.get_good_fiberstatus(zcat)
    good_spec &= zcat['OBJTYPE']=='TGT'    # not included in LSS BGS,LRG,ELG cuts
    zqual['GOOD_SPEC'] = good_spec.copy()  # GOOD_SPEC: true if it is a science spectrum with good hardware status

    for col in ['GOOD_Z_BGS', 'GOOD_Z_LRG', 'GOOD_Z_ELG', 'GOOD_Z_QSO', 'GOOD_Z_LYA']:
        zqual[col] &= zqual['GOOD_SPEC']  # require good hardware quality for GOOD_Z_TRACER

    # Require primary tracer targeting
    if survey in ['main', 'sv1', 'sv2', 'sv3']:
        if survey=='main':
            desi_target_col = 'DESI_TARGET'
        else:
            desi_target_col = survey.upper()+'_DESI_TARGET'

        # The BGS_ANY, LRG+LGE, ELG and QSO target bits are the same in SV1 to main
        is_bgs = (zcat[desi_target_col] & desi_mask.BGS_ANY) != 0
        is_lrg = (zcat[desi_target_col] & (desi_mask.LRG | desi_mask.LGE)) != 0
        is_elg = (zcat[desi_target_col] & desi_mask.ELG) != 0
        is_qso = (zcat[desi_target_col] & desi_mask.QSO) != 0

        # GOOD_Z_TRACER:
        # True if it is a TRACER target and passes TRACER redshift quality cut
        # False if it is not a Tracer target or if it is a TRACER target but fails TRACER redshift quality cut;
        # They apply to the Z column
        zqual['GOOD_Z_BGS'] &= is_bgs
        zqual['GOOD_Z_LRG'] &= is_lrg  # GOOD_Z_LRG includes both LRG and LGE
        zqual['GOOD_Z_ELG'] &= is_elg

        # GOOD_Z_QSO: like GOOD_Z_{BGS,LRG,ELG}, but applies to the Z_QSO column, not the Z column
        # True if it is a QSO target AND passes the QSO redshift quality cut
        zqual['GOOD_Z_QSO'] &= is_qso
        # GOOD_Z_LYA (if available) also applies to the Z_QSO column
        # For GOOD_Z_LYA we do not check for target membership here because it was done in desispec.validredshifts

        # Note that the GOOD_Z_{BGS,LRG,ELG,QSO,LYA} definitions are more restrictive than in desispec.validredshifts
        # as the target membership and GOOD_SPEC requirements are added here

    else:
        for col in ['GOOD_Z_BGS', 'GOOD_Z_LRG', 'GOOD_Z_ELG', 'GOOD_Z_QSO', 'GOOD_Z_LYA']:
            zqual[col] = False

    ######
    # evaluate Z_CONF; proceed from low-confidence to high-confidence

    # default Z_CONF=0 is no confidence
    zqual['Z_CONF'] = np.uint8(0)   # Note: unsigned int because FITS converts signed int8 to bool (!)

    # Z_CONF=1: less confident redshift but maybe ok
    # criteria: the Z_CONF==3 or 2 criteria are not met, but GOOD_SPEC==True & ZWARN==0
    mask = zqual['GOOD_SPEC'] & (zcat['ZWARN']==0)
    zqual['Z_CONF'][mask] = 1

    # Z_CONF=2: placeholder for non-LSS WG supplied quality criteria

    # Z_CONF=3: highly confident redshift
    # criteria: the object must belong to one of the DESI primary extragalactic target classes (BGS, LRG, ELG, QSO)
    # and pass the LSS redshift quality cuts
    mask = zqual['GOOD_Z_BGS'] | zqual['GOOD_Z_LRG'] | zqual['GOOD_Z_ELG'] | zqual['GOOD_Z_QSO'] | zqual['GOOD_Z_LYA']
    zqual['Z_CONF'][mask] = 3

    for col in zqual.colnames:
        if col not in zcat.colnames:
            print(survey, col, 'not in zcat')
        zcat[col] = zqual[col]

    # zcat = hstack([zcat, zqual], join_type='exact')

    # Create "best redshift" columns, choosing between Z and Z_QSO
    z_cols = ['Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'SUBTYPE', 'CHI2', 'DELTACHI2', 'COEFF']
    for col in z_cols:
        zcat[col+'_BEST'] = zcat[col].copy()

    # Use Z_QSO if GOOD_Z_QSO==True and Z_QSO differs by more than 1000 km/s from Z
    c = astropy.constants.c.to('km/s').value
    dv = c*(zcat['Z']-zcat['Z_QSO'])/(1+zcat['Z_QSO'])
    mask = (zcat['GOOD_Z_QSO'] | zcat['GOOD_Z_LYA']) & (np.abs(dv) > 1000)
    zcat['Z_BEST'][mask] = zcat['Z_QSO'][mask].copy()
    for col in z_cols:
        if col!='Z':
            zcat[col+'_BEST'][mask] = zcat[col+'_NEW'][mask].copy()

    cat_stack.append(zcat)

cat_stack = vstack(cat_stack).filled(0)
cat_stack.sort('id')


assert np.all(cat_stack['id'] == cat['id'])
assert np.all(cat_stack['TARGETID'] == cat['TARGETID'])

columns = ['TARGETID', 'SURVEY', 'PROGRAM', 'UNIQPIX', 'Z_BEST', 'Z_CONF', 'ZWARN_BEST', 'SPECTYPE_BEST', 'DELTACHI2_BEST', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'GOOD_SPEC', 'EFFTIME_SPEC', 'ZCAT_PRIMARY', 'DESI_TARGET', 'BGS_TARGET', 'SCND_TARGET', 'GOOD_Z_BGS', 'GOOD_Z_LRG', 'GOOD_Z_ELG', 'GOOD_Z_QSO', 'GOOD_Z_LYA']
cat_stack = cat_stack[columns]

cat_stack.write('/pscratch/sd/r/rongpu/tmp/matterhorn/zcatalog/v2_latest/zall-pix-matterhorn.fits')

