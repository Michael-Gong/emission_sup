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
#E_s=4.1e5
#gg = 100.

eta_list = np.logspace(-5,1,500)
#chi_min  = np.logspace(-19,-7,100)
chi_min  = np.logspace(np.log10(5e-18),np.log10(5e-7),500)
print(chi_min)

f_qed = open('hsokolov.qed', 'w')
f_qed.write('500\t -5\t 1\n')  # python will convert \n to os.linesep
f_rr = open('hsokolov.rr', 'w')
f_rr.write('500\t -5\t 1\n')  # python will convert \n to os.linesep
#f_qed.close()

for i in range(np.size(eta_list)):
    eta = eta_list[i]
#    bin_chi  = np.logspace(np.log10(1e-8*eta),np.log10(0.499999999*eta),2000)
#    grid_chi = np.logspace(np.log10(1e-8*eta),np.log10(0.499999999*eta),2001)
    chi_min_i = chi_min[i]
    bin_chi  = np.logspace(np.log10(chi_min_i),np.log10(0.499999999*eta),2000)
    grid_chi = np.logspace(np.log10(chi_min_i),np.log10(0.499999999*eta),2001)
    chi=bin_chi
    y   = 4*chi/(3*eta*(eta-2*chi))
    F_chi = np.zeros_like(chi)
    for j in range(np.size(chi)):
        result = integrate.quad(lambda x: kv(5./3.,x), y[j], 1e3*y[j])
        F_chi[j] = 4*chi[j]**2/eta**2*y[j]*kv(2./3.,y[j])+(1-2*chi[j]/eta)*y[j]*result[0]
    delta_chi = grid_chi[1:]-grid_chi[:-1]
    print(delta_chi)
    print(np.log10(eta),': ',np.log10(np.sum(F_chi/chi*delta_chi)))    
    f_qed.write(str(np.log10(eta))+'\t'+str(np.log10(np.sum(F_chi/chi*delta_chi)))+'\n') 
  
    bin_chi_c = np.logspace(np.log10(chi_min_i),np.log10(1e2*eta**2),2000)
    grid_chi_c = np.logspace(np.log10(chi_min_i),np.log10(1e2*eta**2),2001)
    chi_c = bin_chi_c
    y_c = 4*chi_c/(3*eta*eta)
    F_chi_c = np.zeros_like(chi_c)
    for j in range(np.size(chi_c)): 
        result = integrate.quad(lambda x: kv(5./3.,x), y_c[j], 1e3*y_c[j])
        F_chi_c[j]=y_c[j]*result[0]
    delta_chi_c = grid_chi_c[1:]-grid_chi_c[:-1]
    print(np.log10(eta),': ',np.log10(np.sum(F_chi_c/chi_c*delta_chi_c)))    
    f_rr.write(str(np.log10(eta))+'\t'+str(np.log10(np.sum(F_chi_c/chi_c*delta_chi_c)))+'\n') 

    print('=================================================')

f_qed.close()
f_rr.close()



