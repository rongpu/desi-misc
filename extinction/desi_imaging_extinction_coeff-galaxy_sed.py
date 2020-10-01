# Based on the code from the LSST-testbed paper, with the deprecated specutils.extinction replaced with dust_extinction

# Calculate the extinction coefficient R_b = A_b/E(B-V)_SFD
# Following the procedure in Schlafly el al. 2011 using synthetic stellar spectrum
# Extrapolate the red end of the synthetic stellar spectrum with blackbody spectrum

from __future__ import division, print_function
import sys, os
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy import units as u
from scipy.interpolate import interp1d
from scipy.integrate import quad

sys.path.append('/Users/rongpu/git/Python/user_modules/')
from user_common import extrap1d

# from specutils.extinction import extinction
from dust_extinction.parameter_averages import F99

r_v = 3.1
ebv = 0.03 # there is a slight dependence of r_b on the E(B-V) value
redshift = 1.

bands = ['decam_g', 'decam_r', 'decam_z', '90prime_g', '90prime_r', 'mosaic_z', 'mosaic_z_corr']
# bands = ['DECam_g', 'DECam_r', 'DECam_z'] # EAZY filter curves
filenames = ['decam.g.am1p4.dat.txt', 'decam.r.am1p4.dat.txt', 'decam.z.am1p4.dat.txt', 'bass-g.txt', 'bass-r.txt', 'kpzd.txt', 'kpzdccdcorr3.txt']

print(bands)

for sed_index in range(1, 8):

    # Galaxy spectrum f_lambda (in flux unit)
    fn = '/Users/rongpu/git/desi-misc/extinction/templates/EAZY_v1.1_lines/eazy_v1.1_sed{}.dat'.format(sed_index)
    # fn = '/Users/rongpu/git/desi-misc/extinction/templates/BR07/default_sed{}.dat'.format(sed_index)
    if os.path.isfile(fn):
        tmp = Table.read(fn, format='ascii')
    else:
        break
    sx, sy = tmp['col1'], tmp['col2']
    # Convert to photon unit
    sy *= sx
    # redshifting
    sx *= (1+redshift)

    # Interpolate stellar spectrum
    s_interp = extrap1d(sx, sy)

    r_b = np.zeros(len(bands))

    for index in range(len(bands)):

        band = bands[index]
        # print(band+' band')

        #---#---#---#---#---#---#---#---#---#---#---#---#---

        # Load filter response curve
        filename = filenames[index]
        tmp = Table.read('/Users/rongpu/git/desi-misc/extinction/filter_curves/'+filename, format='ascii')
        wx = tmp['col1']
        if 'decam' or 'mosaic' in band:
            wy = tmp['col2']
        else:
            wy = tmp['col4']

        # # EAZY filter curves
        # with open('/Users/rongpu/git/desi-misc/extinction/filter_curves/eazy/FILTER.RES.latest', 'r') as f:
        #     while True:
        #         text = f.readline()
        #         if band in text:
        #             nlines = int(text.split()[0])
        #             wx, wy = np.zeros(nlines), np.zeros(nlines)
        #             for ii in range(nlines):
        #                 text = f.readline()
        #                 wx[ii] = float(text.split()[1])
        #                 wy[ii] = float(text.split()[2])
        #             break
        #         if text=='':
        #             print('Error: {} not found!'.format(band))
        #             break

        # Interpolate filter response curve
        w_interp = extrap1d(wx, wy)

        #---#---#---#---#---#---#---#---#---#---#---#---#---

        # Integral in the denominator
        xx = np.linspace(wx.min(), wx.max(), 20000)
        wyy = w_interp(xx)
        syy = s_interp(xx)
        ws = wyy * syy
        ws_interp = extrap1d(xx, ws)
        denominator = quad(ws_interp, xx.min(), xx.max())[0]

        ext = F99(Rv=r_v)

        # extinction at 1um from O'Donnell 94
        # arguments: extinction(wave, a_v, r_v=3.1, model='od94')
        # a_1um_od94 = extinction(np.array([10000])*u.angstrom, r_v*ebv, r_v=r_v, model='od94')[0]
        a_1um_od94 = ebv*1.319
        # extinction at 1um from Fitzpatrick 99
        a_1um_f99 = -2.5*np.log10(ext.extinguish(np.array([10000])*u.angstrom, Ebv=ebv))

        # Integral in the numerator
        a_lambda = -2.5*np.log10(ext.extinguish(xx*u.angstrom, Ebv=ebv))/a_1um_f99
        f = wyy*syy*10.**(-a_lambda*0.78*a_1um_od94/2.5)
        f_interp = extrap1d(xx, f)
        numerator = quad(f_interp, xx.min(), xx.max())[0]

        r_b[index] = -2.5*np.log10(numerator/denominator)/ebv
        # print(r_b[index])

    print('SED {}:'.format(sed_index), r_b)

# Results:

# EAZY_v1.1_lines:
#         decam_g,   decam_r,   decam_z,   90prime_g, 90prime_r, mosaic_z,  mosaic_z_corr
# SED 1: [2.98491469 2.0846769  1.18041201 3.1564938  2.10648058 1.16549885 1.17658546]
# SED 2: [3.24101788 2.16501736 1.20562562 3.30044097 2.17653955 1.18121229 1.19344717]
# SED 3: [3.18457819 2.13810247 1.19775524 3.25524854 2.15346039 1.17713707 1.18913214]
# SED 4: [3.04356625 2.10151823 1.19091625 3.18194708 2.12150032 1.17250961 1.18393599]
# SED 5: [2.89294106 2.08172722 1.18329652 3.10311904 2.10352818 1.16604397 1.17706237]
# SED 6: [3.11888619 2.12006053 1.18230761 3.21037265 2.13922145 1.16246536 1.17472287]
# SED 7: [2.51525564 2.06100682 1.17945711 2.94179305 2.08483511 1.16287941 1.17360909]

# BR07:
#         decam_g,   decam_r,   decam_z,   90prime_g, 90prime_r, mosaic_z,  mosaic_z_corr
# SED 1: [2.38045802 2.03890157 1.17783179 3.04054121 2.08267726 1.16025195 1.1713714 ]
# SED 2: [3.25574971 2.17075813 1.2111619  3.31327209 2.1810943  1.18026935 1.19369327]
# SED 3: [3.13877303 2.11556617 1.19637704 3.24309197 2.13426588 1.17700466 1.18871415]
# SED 4: [2.80936771 2.2720914  1.18467779 3.07872944 2.0986032  1.16697758 1.12881477]
# SED 5: [3.1474041  2.15058524 1.20441921 3.2128008  2.16511195 1.18368875 1.19708353]
