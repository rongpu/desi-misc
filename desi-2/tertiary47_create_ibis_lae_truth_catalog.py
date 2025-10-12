from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio


def read_redrock_catalogs(dirname):
    
    fns = glob.glob(os.path.join(dirname, 'redrock*.fits'))
    # print(len(fns))
    cat = []
    for fn in fns:
        tmp1 = Table(fitsio.read(fn, ext='REDSHIFTS'))
        tmp2 = Table(fitsio.read(fn, ext='FIBERMAP'))
        tmp3 = Table(fitsio.read(fn, ext='TSNR2'))
        assert np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp3['TARGETID'])
        tmp2.remove_column('TARGETID')
        tmp3 = tmp3[['TSNR2_ELG', 'TSNR2_LYA', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG']]
        tmp1['fn'] = os.path.basename(fn)
        cat.append(hstack([tmp1, tmp2, tmp3]))
    cat = vstack(cat)
    # print(len(cat))
    
    return cat


def add_synthg(d, ibis_bands=["M411", "M438", "M464", "M490", "M517"]):
    '''
    Modified from Anand's scripot.
    '''

    cs = np.array([-0.004, 0.120, 0.262, 0.229, 0.359])
    assert len(cs) == len(ibis_bands)

    d["FLUX_SYNTHG"], d["FLUX_IVAR_SYNTHG"], d["FIBERFLUX_SYNTHG"] = 0.0, 0.0, 0.0
    vs = np.zeros(len(d))
    for quant in ["FLUX", "FIBERFLUX"]:
        fkey = "{}_SYNTHG".format(quant)
        for c, band in zip(cs, ibis_bands):
            d[fkey] += c * d["{}_{}".format(quant, band)]
            if quant == "FLUX":
                vs += c**2 / d["FLUX_IVAR_{}".format(band)]
        d["FLUX_IVAR_SYNTHG"] = vs**-1.0
        d["SNR_SYNTHG"] = d["FLUX_SYNTHG"] * np.sqrt(d["FLUX_IVAR_SYNTHG"])

    d.meta["SYNTGCOEF"] = ",".join(cs.astype(str))

    return d


# Load redrock catalog
version = 'laelbg_templates_2.89_3.19_archetype_galaxy_only'
dirname = os.path.join('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary47', version)
cat = read_redrock_catalogs(dirname)
print(len(cat))

mask = cat['COADD_FIBERSTATUS']==0
cat = cat[mask]
print(len(cat))

mask = cat['SCND_TARGET']>0
cat = cat[mask]
print(len(cat))

# Fiber assignment catalog
tt = Table(fitsio.read('/global/cfs/cdirs/desi/survey/fiberassign/special/tertiary/0047/tertiary-targets-0047-assign.fits'))
cat = join(cat, tt[['IBIS', 'IBIS_ROW', 'TARGETID']], keys='TARGETID', join_type='left')
print(len(cat))

# Only keep IBIS-selected targets
mask = cat['IBIS'].copy()
cat = cat[mask]
print(len(cat))

# Confirmed LAEs
lae_tids = [39089837587825585, 39089837587825950, 39089837587826182, 39089837587826547, 39089837587826733, 39089837587826947, 39089837587827165, 39089837587827166, 39089837587827180, 39089837587827434, 39089837587827469, 39089837587827606, 39089837587828052, 39089837587831005, 39089837587831240, 39089837587832640, 39089837587834177, 39089837587834279, 39089837587834451, 39089837587836532, 39089837587836702, 39089837587836776, 39089837587840114, 39089837587841151, 39089837587841597, 39089837587841892, 39089837587843684, 39089837587843780, 39089837587843983, 39089837587844148, 39089837587844511, 39089837587845067, 39089837587845206, 39089837587845256, 39089837587846436, 39089837587850227, 39089837587850803, 39089837587851519, 39089837587851979, 39089837587852083, 39089837587853100, 39089837587854728, 39089837587856510, 39089837587856728, 39089837587856904, 39089837587857046, 39089837587859050, 39089837587859452, 39089837587861292, 39089837587865441, 39089837587865490, 39089837587865856, 39089837587865988, 39089837587866249, 39089837587866834, 39089837587867440, 39089837587868381, 39089837587868635, 39089837587868709, 39089837587870512, 39089837587871317, 39089837587871538, 39089837587871570, 39089837587872590, 39089837587875824, 39089837587877324, 39089837587877425, 39089837587877629, 39089837587880957, 39089837587881267, 39089837587881423, 39089837587881961, 39089837587882264, 39089837587883009, 39089837587884409, 39089837587884433, 39089837587884823, 39089837587886627, 39089837587886758, 39089837587887674, 39089837587888703, 39089837587888900, 39089837587891089, 39089837587891122, 39089837587891186, 39089837587891287, 39089837587892862, 39089837587892997, 39089837587893050, 39089837587893656, 39089837587894583, 39089837587895827, 39089837587895982, 39089837587895999, 39089837587896149, 39089837587896177, 39089837587896889, 39089837587897041, 39089837587897179, 39089837587897412, 39089837587897748, 39089837587898287, 39089837587898382, 39089837587898532, 39089837587900727, 39089837587900911, 39089837587902259, 39089837587902936, 39089837587902996, 39089837587903317, 39089837587903433, 39089837587903468, 39089837587903521, 39089837587903834, 39089837587903836, 39089837587904972, 39089837587905727, 39089837587905885, 39089837587905918, 39089837587907125, 39089837587907130, 39089837587907673, 39089837587908753, 39089837587911660, 39089837587912051, 39089837587824652, 39089837587827367, 39089837587831970, 39089837587839381, 39089837587844279, 39089837587876038, 39089837587882310, 39089837587887716, 39089837587888052, 39089837587897017, 39089837587897893, 39089837587907644, 39089837587911027, 39089837587871407, 39089837587835366, 39089837587906898, 39089837587907193, 39089837587842105, 39089837587855034, 39089837587865676, 39089837587896741, 39089837587897840, 39089837587828057, 39089837587894264, 39089837587890817, 39089837587865038, 39089837587880630, 39089837587881603, 39089837587850895, 39089837587841744, 39089837587882802]
print(len(lae_tids), len(np.unique(lae_tids)))

mask = np.in1d(cat['TARGETID'], lae_tids)
cat['lae'] = mask.copy()
print(np.sum(cat['lae']))

ibis = Table(fitsio.read('/global/cfs/cdirs/desi/survey/fiberassign/special/tertiary/0047/inputcats/ibis4-cosmos-m490-targets-wide.fits'))
print(len(ibis))
ibis['IBIS_ROW'] = np.arange(len(ibis))

# Remove blank columns
cat.remove_columns(['PMRA', 'PMDEC', 'REF_EPOCH', 'RELEASE', 'BRICKNAME', 'BRICKID', 'EBV', 'MASKBITS', 'SERSIC', 'SHAPE_R', 'SHAPE_E1', 'SHAPE_E2', 'REF_ID', 'REF_CAT', 'GAIA_PHOT_G_MEAN_MAG', 'GAIA_PHOT_BP_MEAN_MAG', 'GAIA_PHOT_RP_MEAN_MAG', 'PARALLAX'])
cat.remove_columns(['BRICK_OBJID', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z', 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'PHOTSYS'])
print(len(cat))
cat = join(cat, ibis, keys='IBIS_ROW', join_type='left')
print(len(cat))

cat = add_synthg(cat)  # Add IBIS synthetic g magnitude

with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    cat['m411mag'] = 22.5 - 2.5*np.log10(cat['FLUX_M411']) - 4.290 * cat['EBV']
    cat['m438mag'] = 22.5 - 2.5*np.log10(cat['FLUX_M438']) - 4.099 * cat['EBV']
    cat['m464mag'] = 22.5 - 2.5*np.log10(cat['FLUX_M464']) - 3.877 * cat['EBV']
    cat['m490mag'] = 22.5 - 2.5*np.log10(cat['FLUX_M490']) - 3.634 * cat['EBV']
    cat['m517mag'] = 22.5 - 2.5*np.log10(cat['FLUX_M517']) - 3.389 * cat['EBV']
    cat['m411fibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_M411']) - 4.290 * cat['EBV']
    cat['m438fibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_M438']) - 4.099 * cat['EBV']
    cat['m464fibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_M464']) - 3.877 * cat['EBV']
    cat['m490fibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_M490']) - 3.634 * cat['EBV']
    cat['m517fibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_M517']) - 3.389 * cat['EBV']
    cat['gmag'] = 22.5 - 2.5*np.log10(cat['FORCED_FLUX_G']) - 3.240 * cat['EBV']
    cat['rmag'] = 22.5 - 2.5*np.log10(cat['FORCED_FLUX_R2']) - 2.276 * cat['EBV']
    cat['imag'] = 22.5 - 2.5*np.log10(cat['FORCED_FLUX_I2']) - 1.633 * cat['EBV']
    cat['zmag'] = 22.5 - 2.5*np.log10(cat['FORCED_FLUX_Z']) - 1.263 * cat['EBV']
    cat['ymag'] = 22.5 - 2.5*np.log10(cat['FORCED_FLUX_Y']) - 1.075 * cat['EBV']
    cat['synthgmag'] = 22.5 - 2.5*np.log10(cat['FLUX_SYNTHG']) - 3.240 * cat['EBV']
    cat['synthgfibermag'] = 22.5 - 2.5*np.log10(cat['FIBERFLUX_SYNTHG']) - 3.240 * cat['EBV']

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary47/catalogs/tertiary47_ibis_lae_truth-20250914.fits')
