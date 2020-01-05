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
        'size'   : 28,
       }
font_size = 28
font_size2= 18
def pxpy_to_energy_qed(gamma, weight, eta_i):
    binsize = 50
    en_grid = np.logspace(np.log10(1e-6*eta_i**2),np.log10(0.5*eta_i),50)
    en_bin  = np.logspace(np.log10(1e-6*eta_i**2),np.log10(0.5*eta_i),51)
    en_value = np.zeros_like(en_grid)
    delta_bin = en_bin[1:]-en_bin[:-1] 
    for i in range(binsize):
        en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/delta_bin[i]
    return (en_grid, en_value)

def pxpy_to_energy_rr(gamma, weight, eta_i):
    binsize = 50
    en_grid = np.logspace(np.log10(1e-6*eta_i**2),np.log10(1e1*eta_i**2),50)
    en_bin  = np.logspace(np.log10(1e-6*eta_i**2),np.log10(1e1*eta_i**2),51)
    en_value = np.zeros_like(en_grid)
    delta_bin = en_bin[1:]-en_bin[:-1] 
    for i in range(binsize):
        print(i)
        en_value[i] = sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/delta_bin[i]
    return (en_grid, en_value)

insert1='./qe_15_full/'
px=96.5926 # np.loadtxt(insert1+'px_0000.txt')
py=25.8819 #np.loadtxt(insert1+'py_0000.txt')
gg=(px**2+py**2+1)**0.5
b0=10.0 #np.loadtxt(insert1+'bz_part_0000.txt')
t0=10.472626 #np.loadtxt(insert1+'t_0000.txt')
insert_n='_0000'
photon_px=np.loadtxt(insert1+'photon_px'+insert_n+'.txt')
photon_py=np.loadtxt(insert1+'photon_py'+insert_n+'.txt')
photon = (photon_px**2+photon_py**2)**0.5
chi   = photon*b0/2.0/4.1e5
eta   = gg*b0/4.1e5
dt    = t0/6.28
r0=gg/b0
omega_critic=1.5*gg**3/r0

weight_photon = photon/dt/1e4


ppx,ppy = pxpy_to_energy_rr(chi, weight_photon, eta)
plt.subplot(1,1,1)
#plt.plot(ppx*4*omega_critic/3/eta**2, ppy*3*eta**2/4/omega_critic,  linestyle='-', color='r', marker='o',linewidth=3,markersize=3, label='epoch_qed')
grid_x = ppx*4*omega_critic/3/eta**2
grid_y = ppy*3*eta**2/4/omega_critic
xi =  grid_x*r0/3.0/gg**3
#plt.bar(grid_x, grid_y, width=0.3*grid_x, color='cyan',alpha=0.4,edgecolor='black',linewidth=2,label='Monte-carlo method')
for i in range(np.size(grid_x)):
    result = integrate.quad(lambda x: kv(1.6666666667,x), 2*xi[i], 1e3*2*xi[i])
    grid_y[i] = grid_y[i]*(2*xi[i]*gg*3**0.5/(4*pi**2)*kv(0.66666666666667,xi[i])**2)/result[0]
#plt.plot(grid_x,grid_y,  linestyle='-', color='black', marker='o',linewidth=3,markersize=3, label='epoch_qed')
plt.bar(grid_x, grid_y, width=0.3*grid_x, color='tomato',alpha=0.4,edgecolor='black',linewidth=2,label='Monte-carlo method')


#eta = 0.1
#chi_c = np.logspace(np.log10(1e-5*eta**2),np.log10(1e1*eta**2),2000)
#y_c = 4*chi_c/(3*eta*eta)
#F_chi_c = np.zeros_like(chi_c)
#for i in range(np.size(chi)):
#    result = integrate.quad(lambda x: kv(1.6666666667,x), y_c[i], 1e3*y_c[i])
#    F_chi_c[i]=y_c[i]*result[0]
#F_chi_s = 8*chi**2/3/(3)**0.5/np.pi/eta**4*((1+(1-2*chi/eta)**(-2))*0.921*(y)**(-5.0/3) + 2*(1-2*chi/eta)**(-1)*0.307*(y)**(-5.0/3) + (2*chi/eta)**2*(1-2*chi/eta)**(-2)*0.307*(y)**(-5.0/3)  )
#plt.plot(chi_c**4*omega_critic/3/eta**2, F_chi_c, '--k',  linewidth=3, label=r'$f_{synch}\approx y_c\int_{y_c}^{\infty}K_{5/3}(u)du;\ y_c=\frac{4\chi}{3\eta^2}$')





grid_omega_x = np.logspace(0,np.log10(100*omega_critic),1000)
norm_fac = q0**2/(4*pi*epsilon0*4*pi**2*v0)/(m0*v0**2)*frequency
#data_omega = data_omega*norm_fac
theta_1 = 0
xi =  grid_omega_x*r0/3.0/gg**3*(1.0+gg**2*theta_1**2)**1.5
#y_line = np.linspace(1e-7,5e-4,1000)
#x_line = np.zeros_like(y_line)+omega_critic
theory_line = norm_fac*3*(2*grid_omega_x*r0/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi)**2)
theory_line_1 = np.zeros_like(theory_line)
for i in range(np.size(xi)):
    result = integrate.quad(lambda x: kv(1.6666666667,x), 2*xi[i], 1e3*2*xi[i])
    theory_line_1[i] = theory_line[i]/(2*xi[i]*gg*3**0.5/(4*pi**2)*kv(0.66666666666667,xi[i])**2)*result[0]

#norm_x = matplotlib.colors.Normalize()
plt.plot(grid_omega_x,theory_line,':',color='lime',linewidth=6,label='Synchrotron radiation')
#plt.plot(grid_omega_x,theory_line_1,'--',color='k',linewidth=2,label='S')
#plt.plot(grid_omega_x,theory_line_1,':k',linewidth=2,label='theoretical equation')
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)

plt.grid(which='major',color='k', linestyle=':', linewidth=1)
plt.legend(loc='best',fontsize=font_size2+4,framealpha=0.0)
plt.xscale('log')
plt.yscale('log')
plt.xticks([1e1,1e3,1e5,1e7],fontsize=font_size); plt.yticks(fontsize=font_size);
#plt.xlim(1e-5*eta**2,1e1*eta**2)
plt.ylim(2e-8,6e-4)
plt.xlim(1.5e0,5e7)
#plt.text(1e-4*eta**2,0.6,r'$\eta=0.1$',fontdict=font)

plt.subplots_adjust(top=0.98, bottom=0.15, left=0.15, right=0.98, hspace=0.2, wspace=0.2)
fig = plt.gcf()
fig.set_size_inches(8.4, 8)
#fig.set_size_inches(5, 4.5)
fig.savefig(insert1+'wrap_full_qe15.png',format='png',dpi=160)
plt.close("all")
