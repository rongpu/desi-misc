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


columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'TILEID', 'FIBER', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'PRIORITY', 'OBJTYPE']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG', 'TSNR2_ELG_B', 'TSNR2_BGS_B', 'TSNR2_QSO_B', 'TSNR2_LRG_B']
columns_emline = ['TARGETID', 'OII_FLUX', 'OII_FLUX_IVAR', 'HDELTA_FLUX', 'HDELTA_FLUX_IVAR', 'HGAMMA_FLUX', 'HGAMMA_FLUX_IVAR', 'HBETA_FLUX', 'HBETA_FLUX_IVAR', 'OIII_FLUX', 'OIII_FLUX_IVAR', 'HALPHA_FLUX', 'HALPHA_FLUX_IVAR']

tile_list = [83579, 83580, 83581, 83582]
lastnight_list = [20260116, 20260115, 20260116, 20260116]

###################### Load daily redrock catalogs ######################

print('daily')

fns = glob.glob('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/subsets/*/redrock-*.fits')
fns.sort()
print(len(fns))

cat_stack = []

for fn in fns:
    tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
    tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
    tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
    assert (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID']))
    tmp2.remove_column('TARGETID')
    tmp4.remove_column('TARGETID')
    cat = hstack([tmp1, tmp2, tmp4])
    # cat['LASTNIGHT'] = np.array(os.path.basename(os.path.dirname(fn)), dtype=int)
    cat['coadd_fn'] = fn.replace('redrock-', 'coadd-')

    fn_emline = fn.replace('redrock-', 'emline-')
    emline = Table(fitsio.read(fn_emline, columns=columns_emline))
    emline.remove_column('TARGETID')
    cat = hstack([cat, emline])
    
    cat_stack.append(cat)
    
daily = vstack(cat_stack)

###################### Load custom LAE redrock catalogs ######################

print('LAE NMF')

fns = glob.glob('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/subsets/*/lae_nmf/redrock-*.fits')
fns.sort()
print(len(fns))

cat_stack = []

for fn in fns:
    tmp1 = Table(fitsio.read(fn, ext=1, columns=columns_1))
    tmp2 = Table(fitsio.read(fn, ext=2, columns=columns_2))
    tmp4 = Table(fitsio.read(fn, ext=4, columns=columns_4))
    assert (np.all(tmp1['TARGETID']==tmp2['TARGETID']) and np.all(tmp1['TARGETID']==tmp4['TARGETID']))
    tmp2.remove_column('TARGETID')
    tmp4.remove_column('TARGETID')
    cat = hstack([tmp1, tmp2, tmp4])
    # cat['LASTNIGHT'] = np.array(os.path.basename(os.path.dirname(fn)), dtype=int)
    cat['rr_fn'] = fn

    cat['subset'] = fn.split('/')[-2]

    # fn_emline = fn.replace('redrock-', 'emline-')
    # emline = Table(fitsio.read(fn_emline, columns=columns_emline))
    # emline.remove_column('TARGETID')
    # cat = hstack([cat, emline])
    
    cat_stack.append(cat)
    
cat = vstack(cat_stack)

################################################################################

assert np.all(cat['TARGETID']==daily['TARGETID'])

cat['coadd_fn'] = daily['coadd_fn']

mask = cat['COADD_FIBERSTATUS']==0
print(np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]
daily = daily[mask]

mask = (cat['OBJTYPE']=='TGT') & (cat['PRIORITY']==7000)  # select LAE targets
print(np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]
daily = daily[mask]

cat['EFFTIME_LRG'] = cat['TSNR2_LRG'] * 12.15

################################################################################

# Use daily emission line S/N to identify low-z galaxies
daily['oii_sn'] = daily['OII_FLUX'] * np.sqrt(daily['OII_FLUX_IVAR'])
daily['halpha_sn'] = daily['HALPHA_FLUX'] * np.sqrt(daily['HALPHA_FLUX_IVAR'])
daily['hbeta_sn'] = daily['HBETA_FLUX'] * np.sqrt(daily['HBETA_FLUX_IVAR'])
daily['oiii_sn'] = daily['OIII_FLUX'] * np.sqrt(daily['OIII_FLUX_IVAR'])

mask_star = (daily['SPECTYPE']=='STAR') & (daily['ZWARN']==0)

mask_z = daily['Z']<0.49
mask_oii = daily['oii_sn'] > 3.
mask_halpha = daily['halpha_sn'] > 3.
mask_oiii = daily['oiii_sn'] > 4.
mask_halpha_1 = daily['halpha_sn'] > 5.
mask_oiii_1 = daily['oiii_sn'] > 5.
mask = mask_z & (mask_oii & (mask_halpha | mask_oiii) | (mask_halpha_1 & mask_oiii_1))
print(np.sum(mask))
mask1 = mask.copy()

mask_z = (daily['Z']>0.4) & (daily['Z']<0.95)
mask_oii = daily['oii_sn'] > 3.
mask_oiii = daily['oiii_sn'] > 3.
mask = mask_z & (mask_oii & mask_oiii)
print(np.sum(mask))
mask2 = mask.copy()

mask_lowz = mask_star | mask1 | mask2
print('mask_lowz', np.sum(mask_lowz))

mask_z = (daily['Z']>0.8)  & (daily['Z']<1.6)
mask_elg = mask_z & (daily['OII_FLUX']>0) & (daily['OII_FLUX_IVAR']>0)
mask_elg &= np.log10(daily['OII_FLUX'] * np.sqrt(daily['OII_FLUX_IVAR'])) > 0.9 - 0.2 * np.log10(daily['DELTACHI2'])
mask_elg &= cat['DELTACHI2']<20  # override the ELG solution if the LAE fit is confident; almost all the overrode redshifts are fitting sky residuals
print('mask_elg', np.sum(mask_elg))

mask_qso = daily['ZWARN']==0
mask_qso &= daily['DELTACHI2']>100
z_diff_threshold = 0.01  # 3000 km/s
mask_qso &= np.abs((daily['Z'] - cat['Z'])/(1 + cat['Z'])) > z_diff_threshold
mask_qso &= daily['SPECTYPE']=='QSO'
# mask_qso &= daily['Z']<2.25  # only remove low-z QSOs; not implemented due to negative flux artifacts being fitted by QSO templates

mask_contam = (mask_lowz | mask_elg | mask_qso)
print('mask_contam', np.sum(mask_contam))


# DELTACHI2 threshold that increases with redshift
def dchi2_vs_z(z):
    min_dchi2 = np.zeros(len(z))
    mask = z<=2.8
    min_dchi2[mask] = 9
    mask = z>2.8
    min_dchi2[mask] = 10**(np.log10(9) + (np.log10(18)-np.log10(9))/0.6 * (z[mask]-2.8))
    return min_dchi2

mask_lae = (cat['Z']>2.2) & (cat['ZWARN']==0) & (cat['DELTACHI2']>dchi2_vs_z(cat['Z'])) & (~mask_contam)
print('LAEs', np.sum(mask_lae), np.sum(mask_lae)/len(mask_lae))

# print('Remove known low-z interlopers')
# interlopers = np.loadtxt('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/catalogs/known_lowz_interlopers.txt', dtype='int')
# mask_lae &= ~np.in1d(cat['TARGETID'], interlopers)
# print('LAEs', np.sum(mask_lae), np.sum(mask_lae)/len(mask_lae))

cat['low_z'] = mask_contam.copy()
cat['lae'] = mask_lae.copy()

cat['Z_daily'] = daily['Z']
cat['ZWARN_daily'] = daily['ZWARN']
cat['DELTACHI2_daily'] = daily['DELTACHI2']
cat['SPECTYPE_daily'] = daily['SPECTYPE']

cat.write('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/catalogs/tertiary49_ibis_lae_subsets-20260127.fits', overwrite=False)

