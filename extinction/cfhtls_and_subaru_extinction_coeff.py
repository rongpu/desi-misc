# Based on the code from the LSST-testbed paper, with the deprecated specutils.extinction replaced with dust_extinction

# Calculate the extinction coefficient R_b = A_b/E(B-V)_SFD
# Following the procedure in Schlafly el al. 2011 using synthetic stellar spectrum
# Extrapolate the red end of the synthetic stellar spectrum with blackbody spectrum

from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
from astropy import units as u
from scipy.interpolate import interp1d
from scipy.integrate import quad
import sys
sys.path.append('/Users/rongpu/git/Python/user_modules/')
from user_common import extrap1d

# from specutils.extinction import extinction
from dust_extinction.parameter_averages import F99

r_v = 3.1

bands = ['u', 'g', 'r', 'i', 'i2', 'z', 'y']
filenames = ['uMega.fil', 'gMega.fil', 'rMega.fil', 'iMega.fil', 'i2Mega.fil', 'zMega.fil', 'suprime_y.txt']

# Stellar spectrum
# In flux unit
sy = np.loadtxt('/Users/rongpu/git/LSST-Testbed/Others/extinction/Munari 7000K spectrum/T07000G45M10V000K1AODNVD01F.ASC')
sx = np.linspace(2500.5, 10499.5, 8000)
# Convert to photon unit
sy1 = sy*sx

sx_extrap = np.linspace(10500.5, 13000.5, 2500)

# # Rayleigh-Jeans extrapolation
# rj_law = sx_extrap**(-3) * (sy1[-10]/sx[-10]**(-3))

# Planck's law * normalization
h = 6.62607004*10**(-34)
c = 3.*10**8
k = 1.38064852*10**(-23)
t = 7000.
sx_extrap = np.linspace(10500.5, 13000.5, 2500)
bb = sx_extrap**(-4)/(np.exp(h*c/(sx_extrap*10**(-10)*k*t))-1) * (sy1[-10]/(sx[-10]**(-4)/(np.exp(h*c/(sx[-10]*10**(-10)*k*t))-1)))

# Combine two arrays
sx_combined = np.concatenate((sx, sx_extrap))
# sy_combined = np.concatenate((sy1, rj_law))
sy_combined = np.concatenate((sy1, bb))

# plt.plot(sx_combined, sy_combined)
# plt.show()

# Interpolate stellar spectrum
# s_interp = extrap1d(sx, sy1)
s_interp = extrap1d(sx_combined, sy_combined)

r_b = np.zeros(len(bands))

for index in range(len(bands)):

    band = bands[index]
    filename = filenames[index]
    print(band+' band')

    #---#---#---#---#---#---#---#---#---#---#---#---#---

    # Load filter response curve

    # CFHTLS total filter curves are downloaded from here:
    # http://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/en/megapipe/docs/filt.html
    # Subaru total filter curve is calculated here:
    # Photo-z/filter_curves/suprime_y.py
    data = np.transpose(np.loadtxt('/Users/rongpu/git/LSST-Testbed/Others/extinction/total filter response/'+filename))
    wx = data[0]
    wy = data[1]

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

    ebv = 0.0001
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

print(bands)
print(r_b)

