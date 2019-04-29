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
        'size'   : 20,
       }
font_size = 20

gg=np.loadtxt('./bench_data/py_0000.txt')
b0=np.loadtxt('./bench_data/bz_part_0000.txt')
t0=np.loadtxt('./bench_data/t_0000.txt')
gg = gg[0]
b0 = b0[1]
t0 = t0[-1]
insert1='./bench_data/'
insert_n='_0000'
photon_px=np.loadtxt(insert1+'photon_px'+insert_n+'.txt')
photon_py=np.loadtxt(insert1+'photon_py'+insert_n+'.txt')
photon = (photon_px**2+photon_py**2)**0.5
chi   = photon*b0/2.0/4.1e5
gamma = gg
eta   = gg*b0/4.1e5
dt    = t0/6.28*3.333e-15
num_bins = 1000
grid_a = np.logspace(np.log10(1e-6),np.log10(0.5), num_bins+1)
weight_photon = gamma/eta*h_planck*chi/1.7/(m0*v0**2)*137/dt/(0.5/num_bins)
#weight_photon = gamma
# the histogram of the data
n, bins, patches = plt.hist(chi, bins=np.linspace(1e-4,0.5, num_bins),weights=weight_photon, range=(1e-4,0.5), normed=0, facecolor='red', alpha=0.6,label='epoch_qed')

gg=np.loadtxt('./bench_data_qed/py_0000.txt')
b0=np.loadtxt('./bench_data_qed/bz_part_0000.txt')
t0=np.loadtxt('./bench_data_qed/t_0000.txt')
gg = gg[0]
b0 = b0[1]
t0 = t0[-1]
insert1='./bench_data_qed/'
insert_n='_0000'
photon_px=np.loadtxt(insert1+'photon_px'+insert_n+'.txt')
photon_py=np.loadtxt(insert1+'photon_py'+insert_n+'.txt')
photon = (photon_px**2+photon_py**2)**0.5
chi   = photon*b0/2.0/4.1e5
gamma = gg
eta   = gg*b0/4.1e5
dt    = t0/6.28*3.333e-15
num_bins = 1000
grid_a = np.logspace(np.log10(1e-6),np.log10(0.5), num_bins+1)
weight_photon = gamma/eta*h_planck*chi/1.7/(m0*v0**2)*137/dt/(0.5/num_bins)
#weight_photon = gamma
# the histogram of the data
n, bins, patches = plt.hist(chi, bins=np.linspace(1e-4,0.5, num_bins),weights=weight_photon, range=(1e-4,0.5), normed=0, facecolor='blue', alpha=0.6,label='my_qed')

gg=np.loadtxt('./bench_data_rr/py_0000.txt')
b0=np.loadtxt('./bench_data_rr/bz_part_0000.txt')
t0=np.loadtxt('./bench_data_rr/t_0000.txt')
gg = gg[0]
b0 = b0[1]
t0 = t0[-1]
insert1='./bench_data_rr/'
insert_n='_0000'
photon_px=np.loadtxt(insert1+'photon_px'+insert_n+'.txt')
photon_py=np.loadtxt(insert1+'photon_py'+insert_n+'.txt')
photon = (photon_px**2+photon_py**2)**0.5
chi   = photon*b0/2.0/4.1e5
gamma = gg
eta   = gg*b0/4.1e5
dt    = t0/6.28*3.333e-15
num_bins = 1000
grid_a = np.logspace(np.log10(1e-4),np.log10(1e1), num_bins+1)
weight_photon = gamma/eta*h_planck*chi/1.7/(m0*v0**2)*137/dt/(1e1/num_bins)
#weight_photon = gamma
# the histogram of the data
n, bins, patches = plt.hist(chi, bins=np.linspace(1e-4,1e1, num_bins),weights=weight_photon, range=(1e-4,1e1), normed=0, facecolor='green', alpha=0.6,label='my_classical')

b0 = 200.0
g0 = 2050.0
eta = g0*b0/4.1e5
#eta = 0.1
chi = np.logspace(np.log10(1e-5*eta**2),np.log10(0.49*eta**2),2000)
chi_c = np.logspace(np.log10(1e-5*eta**2),np.log10(1e2*eta**2),2000)
y   = 4*chi/(3*eta*(eta-2*chi))
y_c = 4*chi_c/(3*eta*eta)
F_chi = np.zeros_like(chi)
F_chi_c = np.zeros_like(chi_c)
for i in range(np.size(chi)):
    result = integrate.quad(lambda x: kv(1.6666666667,x), y[i], 1e3*y[i])
    F_chi[i] = 4*chi[i]**2/eta**2*y[i]*kv(0.6666666667,y[i])+(1-2*chi[i]/eta)*y[i]*result[0]
#    print(result)
    result = integrate.quad(lambda x: kv(1.6666666667,x), y_c[i], 1e3*y_c[i])
    F_chi_c[i]=y_c[i]*result[0]
#F_chi_s = 8*chi**2/3/(3)**0.5/np.pi/eta**4*((1+(1-2*chi/eta)**(-2))*0.921*(y)**(-5.0/3) + 2*(1-2*chi/eta)**(-1)*0.307*(y)**(-5.0/3) + (2*chi/eta)**2*(1-2*chi/eta)**(-2)*0.307*(y)**(-5.0/3)  )
plt.plot(chi, F_chi,   '-k',   linewidth=3, label=r'$F(\eta,\chi)\approx\frac{4\chi^2}{\eta^2}yK_{2/3}(y)+(1-\frac{2\chi}{\eta})y\int_{y}^{\infty}K_{5/3}(t)dt;\ y=\frac{4\chi}{3\eta(\eta-2\chi)}$')
plt.plot(chi_c, F_chi_c, '--k',  linewidth=3, label=r'$f_{synch}\approx y_c\int_{y_c}^{\infty}K_{5/3}(u)du;\ y_c=\frac{4\chi}{3\eta^2}$')

plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.legend(loc='best',fontsize=12,framealpha=0.0)
plt.xticks(fontsize=font_size); plt.yticks(fontsize=font_size);
plt.xlabel(r'$\chi$',fontdict=font)
plt.ylabel(r'$F(\eta,\chi)$',fontdict=font)
plt.xscale('log')
plt.xlim(1e-4,1e1)
#plt.ylim(0,1)
plt.text(1e-3,0.6,r'$\eta=1.0$',fontdict=font)

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)

fig = plt.gcf()
fig.set_size_inches(10, 8.5)
#fig.set_size_inches(5, 4.5)
fig.savefig('./plot_bench_eta=1.png',format='png',dpi=160)
#plt.close("all")
