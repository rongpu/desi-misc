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

sys.path.append(os.path.expanduser('~/git/desi-examples/misc/plot_spectrum'))
from desi_plot_spectrum import plot_spectrum

lyman_alpha_z = {'M411': [2.3, 2.51], 'M438': [2.51, 2.72], 'M464': [2.72, 2.93], 'M490': [2.93, 3.15], 'M517': [3.15, 3.36]}

columns_1 = ['TARGETID', 'CHI2', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'DELTACHI2']
columns_2 = ['TARGETID', 'TILEID', 'FIBER', 'PETAL_LOC', 'DEVICE_LOC', 'LOCATION', 'COADD_FIBERSTATUS', 'TARGET_RA', 'TARGET_DEC', 'MORPHTYPE', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'PARALLAX', 'EBV', 'FLUX_W1', 'FLUX_W2', 'FIBERFLUX_Z', 'MASKBITS', 'PHOTSYS', 'DESI_TARGET', 'BGS_TARGET', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'PRIORITY', 'OBJTYPE']
columns_4 = ['TARGETID', 'TSNR2_ELG', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG', 'TSNR2_ELG_B', 'TSNR2_BGS_B', 'TSNR2_QSO_B', 'TSNR2_LRG_B']
columns_emline = ['TARGETID', 'OII_FLUX', 'OII_FLUX_IVAR', 'HDELTA_FLUX', 'HDELTA_FLUX_IVAR', 'HGAMMA_FLUX', 'HGAMMA_FLUX_IVAR', 'HBETA_FLUX', 'HBETA_FLUX_IVAR', 'OIII_FLUX', 'OIII_FLUX_IVAR', 'HALPHA_FLUX', 'HALPHA_FLUX_IVAR']

tile_list = [83579, 83580, 83581, 83582]
lastnight_list = [20260116, 20260115, 20260116, 20260116]

fns = []
for tileid, lastnight in zip(tile_list, lastnight_list):
    print(tileid, lastnight)
    fns += glob.glob('/global/cfs/cdirs/desi/spectro/redux/daily/tiles/cumulative/{}/{}/redrock-*thru{}.fits'.format(tileid, lastnight, lastnight))
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
    cat['fn'] = fn

    fn_emline = fn.replace('redrock-', 'emline-')
    emline = Table(fitsio.read(fn_emline, columns=columns_emline))
    emline.remove_column('TARGETID')
    cat = hstack([cat, emline])
    
    cat_stack.append(cat)
    
cat = vstack(cat_stack)

print(len(cat), len(np.unique(cat['TARGETID'])))
cat.sort('TSNR2_LRG', reverse=True)
_, idx_keep = np.unique(cat['TARGETID'], return_index=True)
cat = cat[idx_keep]
print(len(cat), len(np.unique(cat['TARGETID'])))

daily = cat.copy()

fns = []
for tileid, lastnight in zip(tile_list, lastnight_list):
    print(tileid, lastnight)
    fns += glob.glob('/global/cfs/cdirs/desicollab/users/rongpu/data/desi2/tertiary49/rr_lae_2.25_3.4_nmf/redrock-*-{}-thru{}.fits'.format(tileid, lastnight))
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
    cat['fn'] = fn

    # fn_emline = fn.replace('redrock-', 'emline-')
    # emline = Table(fitsio.read(fn_emline, columns=columns_emline))
    # emline.remove_column('TARGETID')
    # cat = hstack([cat, emline])
    
    cat_stack.append(cat)
    
cat = vstack(cat_stack)

print(len(cat), len(np.unique(cat['TARGETID'])))
cat.sort('TSNR2_LRG', reverse=True)
_, idx_keep = np.unique(cat['TARGETID'], return_index=True)
cat = cat[idx_keep]
print(len(cat), len(np.unique(cat['TARGETID'])))

assert np.all(cat['TARGETID']==daily['TARGETID'])

mask = cat['COADD_FIBERSTATUS']==0
print(np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]
daily = daily[mask]

mask = (cat['OBJTYPE']=='TGT') & (cat['PRIORITY']==7000)  # select LAE targets
print(np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]
daily = daily[mask]

cat['EFFTIME_LRG'] = cat['TSNR2_LRG'] * 12.15
plt.hist(cat['EFFTIME_LRG'], 100)
plt.axvline(3000, lw=1, ls='--', color='r')
plt.show()

mask = cat['EFFTIME_LRG']>3000
print(np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]
daily = daily[mask]

plt.hist(cat['TSNR2_ELG_B'], 100)
plt.axvline(0.4, lw=1, ls='--', color='r')
plt.show()

mask = cat['TSNR2_ELG_B']>0.4
print(np.sum(mask), np.sum(mask)/len(mask))
cat = cat[mask]
daily = daily[mask]

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

mask = (cat['Z']>2.2) & (cat['ZWARN']==0) & (~mask_contam)
print(np.sum(mask), np.sum(mask)/len(mask))

def dchi2_vs_z(z):
    min_dchi2 = np.zeros(len(z))
    mask = z<=2.8
    min_dchi2[mask] = 9
    mask = z>2.8
    min_dchi2[mask] = 10**(np.log10(9) + (np.log10(18)-np.log10(9))/0.6 * (z[mask]-2.8))
    return min_dchi2

mask_lae = (cat['Z']>2.2) & (cat['ZWARN']==0) & (cat['DELTACHI2']>dchi2_vs_z(cat['Z'])) & (~mask_contam)
print(np.sum(mask_lae), np.sum(mask_lae)/len(mask_lae))

plt.figure(figsize=(7, 7))
plt.plot(cat['Z'][mask_lae], cat['DELTACHI2'][mask_lae], '.', ms=1.5)
plt.xlim(2.2, 3.45)
plt.yscale('log')
# plt.ylim(9, 1e5)
plt.xlabel('Z')
plt.grid(alpha=0.5)
plt.show()

from PIL import Image

def create_vi_png(index):
    tid = cat['TARGETID'][index]
    coadd_fn = daily['fn'][index].replace('redrock-', 'coadd-')
    redrock_fn = cat['fn'][index]
    z = cat['Z'][index]
    
    image_files = []

    plot_path = os.path.join(tmp_dir, str(tid)+'_1.png')
    image_files.append(plot_path)
    plot_spectrum(coadd_fn, tid, redrock_fn=redrock_fn, use_targetid=True, show=False, ylim=1.5, show_model=True, save_path=plot_path)
    plt.close()
    
    plot_path = os.path.join(tmp_dir, str(tid)+'_2.png')
    image_files.append(plot_path)
    plot_spectrum(coadd_fn, tid, redrock_fn=redrock_fn, use_targetid=True, show=False, ylim='minmax', show_model=True, xlim=[(1+z)*1215.67-120, (1+z)*1215.67+120], gauss_smooth=None, figsize=(11, 5), save_path=plot_path, show_text=False)
    plt.close()

    z_oii = (1+z) * 1215.67 / 3727.424 - 1
    plot_path = os.path.join(tmp_dir, str(tid)+'_3.png')
    image_files.append(plot_path)
    plot_spectrum(coadd_fn, tid, z=z_oii, use_targetid=True, show=False, ylim=1.5, show_model=False, save_path=plot_path)
    plt.close()
    
    images = [Image.open(f) for f in image_files]
    heights = [img.height for img in images]
    widths = [img.width for img in images]
    total_height = sum(heights)
    max_width = max(widths)

    stitched = Image.new('RGB', (max_width, total_height), color='white')
    y_offset = 0
    for img in images:
        # Convert to RGB if needed (in case of RGBA or other modes)
        if img.mode != 'RGB':
            img = img.convert('RGB')
        stitched.paste(img, (0, y_offset))
        y_offset += img.height
    
    plot_path = os.path.join(plot_dir, str(tid)+'.png')
    stitched.save(plot_path)
    
    for fn in image_files:
        os.remove(fn)
    
    return plot_path


# tmp_dir = '/pscratch/sd/r/rongpu/tmp/plot_spectrum'
# plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/plots/desi2/tertiary49/lae_candidates_20260124'

# mask = mask_lae.copy()
# idx = np.where(mask)[0]
# print(len(idx))

# print('Plotting')
# from multiprocessing import Pool
# n_process = 128
# with Pool(processes=n_process) as pool:
#     res = pool.map(create_vi_png, idx)


# ########################################################################################################################################################

# tmp_dir = '/pscratch/sd/r/rongpu/tmp/plot_spectrum'
# plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/plots/desi2/tertiary49/lae_candidates_20260124_low_dchi2'

# def dchi2_vs_z1(z):
#     min_dchi2 = np.zeros(len(z))
#     mask = z<=2.8
#     min_dchi2[mask] = 12
#     mask = z>2.8
#     min_dchi2[mask] = 10**(np.log10(12) + (np.log10(20)-np.log10(12))/0.6 * (z[mask]-2.8))
#     return min_dchi2

# mask_lae = (cat['Z']>2.2) & (cat['ZWARN']==0) & (cat['DELTACHI2']>dchi2_vs_z(cat['Z'])) & (~mask_contam)
# mask_lae1 = (cat['Z']>2.2) & (cat['ZWARN']==0) & (cat['DELTACHI2']>dchi2_vs_z1(cat['Z'])) & (~mask_contam)
# mask = mask_lae & (~mask_lae1)
# print(np.sum(mask_lae), np.sum(mask_lae)/len(mask_lae))

# idx = np.where(mask)[0]
# print(len(idx))

# print('Plotting')
# from multiprocessing import Pool
# n_process = 128
# with Pool(processes=n_process) as pool:
#     res = pool.map(create_vi_png, idx)


########################################################################################################################################################

def create_vi_png_for_non_lae(index):
    tid = cat['TARGETID'][index]
    coadd_fn = daily['fn'][index].replace('redrock-', 'coadd-')
    redrock_fn = daily['fn'][index]
    plot_path = os.path.join(plot_dir, str(tid)+'.png')
    plot_spectrum(coadd_fn, tid, redrock_fn=redrock_fn, use_targetid=True, show=False, ylim=1.5, show_model=True, save_path=plot_path)
    plt.close()


###################### Confident low-z interlopers ######################
plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/plots/desi2/tertiary49/lae_candidates_20260124_lowz_interloper_examples'

mask = mask_contam.copy()
idx = np.where(mask)[0]
print(len(idx))

n_plot = 128
if len(idx)>n_plot:
    np.random.seed(666)
    idx = np.random.choice(idx, replace=False, size=n_plot)

print('Plotting')
from multiprocessing import Pool
n_process = 128
with Pool(processes=n_process) as pool:
    res = pool.map(create_vi_png_for_non_lae, idx)

###################### don't know what they are ######################
plot_dir = '/global/cfs/cdirs/desicollab/users/rongpu/plots/desi2/tertiary49/lae_candidates_20260124_unknown_examples'

mask = (~mask_contam) & (~mask_lae)
idx = np.where(mask)[0]
print(len(idx))

n_plot = 128
if len(idx)>n_plot:
    np.random.seed(666)
    idx = np.random.choice(idx, replace=False, size=n_plot)

print('Plotting')
from multiprocessing import Pool
n_process = 128
with Pool(processes=n_process) as pool:
    res = pool.map(create_vi_png_for_non_lae, idx)
