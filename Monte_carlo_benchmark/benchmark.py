import sdf
import matplotlib
matplotlib.use('agg')
#%matplotlib inline
import matplotlib.pyplot as plt
import numpy as np
#from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import scipy.integrate as integrate
import scipy.special as special
from scipy.special import kv
  
######## Constant defined here ########
pi        =     3.1415926535897932384626
pi2d      =     180./pi
q0        =     1.602176565e-19 # C
m0        =     9.10938291e-31  # kg
v0        =     2.99792458e8    # m/s^2
kb        =     1.3806488e-23   # J/K
mu0       =     4.0e-7*pi       # N/A^2
epsilon0  =     8.8541878176203899e-12 # F/m
h_planck  =     6.62606957e-34  # J s
wavelength=     1.0e-6
frequency =     v0*2*pi/wavelength

exunit    =     m0*v0*frequency/q0
bxunit    =     m0*frequency/q0
denunit    =     frequency**2*epsilon0*m0/q0**2
print('electric field unit: '+str(exunit))
print('magnetic field unit: '+str(bxunit))
print('density unit nc: '+str(denunit))

font = {'family' : 'monospace',
        'style'  : 'normal',
        'color'  : 'black',
            'weight' : 'normal',
        'size'   : 25,
       }
font_size = 25

#weight_photon = gamma
# the histogram of the data

plt.subplot(1,2,1)
eta = 0.01
chi = np.logspace(np.log10(1e-5*eta),np.log10(0.499999999*eta),200)
y   = 4*chi/(3*eta*(eta-2*chi))
y_c = 4*chi/(3*eta*eta)
F_chi = np.zeros_like(chi)
F_chi_c = np.zeros_like(chi)
for i in range(np.size(chi)):
    result = integrate.quad(lambda x: kv(5./3.,x), y[i], 1e2*y[i])
    F_chi[i] = 4*chi[i]**2/eta**2*y[i]*kv(2./3.,y[i])+(1-2*chi[i]/eta)*y[i]*result[0]
    result = integrate.quad(lambda x: kv(5./3.,x), y_c[i], 1e2*y_c[i])
    F_chi_c[i]=y_c[i]*result[0]
plt.plot(chi, F_chi,   '-k',   linewidth=3, label='$\eta=0.01$')
plt.plot(chi, F_chi_c, '--k',  linewidth=3)
eta = 0.1
chi = np.logspace(np.log10(1e-6*eta),np.log10(0.499999999*eta),200)
y   = 4*chi/(3*eta*(eta-2*chi))
y_c = 4*chi/(3*eta*eta)
F_chi = np.zeros_like(chi)
F_chi_c = np.zeros_like(chi)
for i in range(np.size(chi)):
    result = integrate.quad(lambda x: kv(5./3.,x), y[i], 1e2*y[i])
    F_chi[i] = 4*chi[i]**2/eta**2*y[i]*kv(2./3.,y[i])+(1-2*chi[i]/eta)*y[i]*result[0]
    result = integrate.quad(lambda x: kv(5./3.,x), y_c[i], 1e2*y_c[i])
    F_chi_c[i]=y_c[i]*result[0]
plt.plot(chi, F_chi,   '-b',   linewidth=3, label='$\eta=0.1$')
plt.plot(chi, F_chi_c, '--b',  linewidth=3)
eta = 1.
chi = np.logspace(np.log10(1e-7*eta),np.log10(0.499999999*eta),200)
y   = 4*chi/(3*eta*(eta-2*chi))
y_c = 4*chi/(3*eta*eta)
F_chi = np.zeros_like(chi)
F_chi_c = np.zeros_like(chi)
for i in range(np.size(chi)):
    result = integrate.quad(lambda x: kv(5./3.,x), y[i], 1e2*y[i])
    F_chi[i] = 4*chi[i]**2/eta**2*y[i]*kv(2./3.,y[i])+(1-2*chi[i]/eta)*y[i]*result[0]
    result = integrate.quad(lambda x: kv(5./3.,x), y_c[i], 1e2*y_c[i])
    F_chi_c[i]=y_c[i]*result[0]
plt.plot(chi, F_chi,   '-r',   linewidth=3, label='$\eta=1$')
plt.plot(chi, F_chi_c, '--r',  linewidth=3)

plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.legend(loc='best',fontsize=20,framealpha=0.2)
plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
plt.xlabel(r'$\chi$',fontdict=font)
plt.ylabel(r'$F(\eta,\chi)$',fontdict=font)
plt.xscale('log')
#plt.yscale('log')
#plt.xlim(1e-4,1e1)
plt.ylim(0,1)
#plt.text(1,0.6,r'$\eta=0.1$',fontdict=font)
# Tweak spacing to prevent clipping of ylabel



plt.subplot(1,2,2)
eta = 0.01
chi = np.logspace(np.log10(1e-5*eta),np.log10(0.499999999*eta),200)
y   = 4*chi/(3*eta*(eta-2*chi))
y_c = 4*chi/(3*eta*eta)
F_chi = np.zeros_like(chi)
F_chi_c = np.zeros_like(chi)
for i in range(np.size(chi)):
    result = integrate.quad(lambda x: kv(5./3.,x), y[i], 1e2*y[i])
    F_chi[i] = 4*chi[i]**2/eta**2*y[i]*kv(2./3.,y[i])+(1-2*chi[i]/eta)*y[i]*result[0]
    result = integrate.quad(lambda x: kv(5./3.,x), y_c[i], 1e2*y_c[i])
    F_chi_c[i]=y_c[i]*result[0]
plt.plot(chi, F_chi/chi,   '-k',   linewidth=3, label='$\eta=0.01$')
plt.plot(chi, F_chi_c/chi, '--k',  linewidth=3)
eta = 0.1
chi = np.logspace(np.log10(1e-6*eta),np.log10(0.499999999*eta),200)
y   = 4*chi/(3*eta*(eta-2*chi))
y_c = 4*chi/(3*eta*eta)
F_chi = np.zeros_like(chi)
F_chi_c = np.zeros_like(chi)
for i in range(np.size(chi)):
    result = integrate.quad(lambda x: kv(5./3.,x), y[i], 1e2*y[i])
    F_chi[i] = 4*chi[i]**2/eta**2*y[i]*kv(2./3.,y[i])+(1-2*chi[i]/eta)*y[i]*result[0]
    result = integrate.quad(lambda x: kv(5./3.,x), y_c[i], 1e2*y_c[i])
    F_chi_c[i]=y_c[i]*result[0]
plt.plot(chi, F_chi/chi,   '-b',   linewidth=3, label='$\eta=0.1$')
plt.plot(chi, F_chi_c/chi, '--b',  linewidth=3)
eta = 1.
chi = np.logspace(np.log10(1e-7*eta),np.log10(0.499999999*eta),200)
y   = 4*chi/(3*eta*(eta-2*chi))
y_c = 4*chi/(3*eta*eta)
F_chi = np.zeros_like(chi)
F_chi_c = np.zeros_like(chi)
for i in range(np.size(chi)):
    result = integrate.quad(lambda x: kv(5./3.,x), y[i], 1e2*y[i])
    F_chi[i] = 4*chi[i]**2/eta**2*y[i]*kv(2./3.,y[i])+(1-2*chi[i]/eta)*y[i]*result[0]
    result = integrate.quad(lambda x: kv(5./3.,x), y_c[i], 1e2*y_c[i])
    F_chi_c[i]=y_c[i]*result[0]
plt.plot(chi, F_chi/chi,   '-r',   linewidth=3, label='$\eta=1$')
plt.plot(chi, F_chi_c/chi, '--r',  linewidth=3)

plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.legend(loc='best',fontsize=20,framealpha=0.2)
plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
plt.xlabel(r'$\chi$',fontdict=font)
plt.ylabel(r'$F(\eta,\chi)/\chi$',fontdict=font)
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-7,1e1)
plt.ylim(1e-1,1e7)
#plt.text(1,0.6,r'$\eta=0.1$',fontdict=font)
# Tweak spacing to prevent clipping of ylabel


plt.subplots_adjust(left=0.12, bottom=0.12, right=0.95, top=0.95,
                wspace=0.3, hspace=0.2)

fig = plt.gcf()
fig.set_size_inches(18, 8.5)
#fig.set_size_inches(5, 4.5)
fig.savefig('./bench_eta=01.png',format='png',dpi=160)
#plt.close("all")
