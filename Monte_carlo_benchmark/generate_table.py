import sdf
import matplotlib
#matplotlib.use('agg')
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


hso_eta = np.linspace(-5,1,100)
hso_dy  = np.zeros_like(hso_eta)
hso_dy_c= np.zeros_like(hso_eta)

hso_chi_min = np.log10(5)+np.linspace(-18,-7,100)

tau_c = 1.05e-34/(0.51*1.6e-13)
for i in range(np.size(hso_eta)):
    log_eta = hso_eta[i]
    log_chi = np.linspace(hso_chi_min[i],np.log10(0.49999)+log_eta,100)   
    d_log_chi = log_chi[-1]-log_chi[-2] 
#    chi_grid = np.linspace(hso_chi_min[i],0.49999*eta,100)    
#    chi = np.linspace(1e-4*eta,0.49*eta,100)    
    eta = 10**log_eta
    chi = 10**log_chi
    y   = 4*chi/(3*eta*(eta-2*chi))
    y_c = 4*chi/(3*eta*eta)
    F_chi   = np.zeros_like(chi)
    F_chi_c = np.zeros_like(chi)
    for j in range(np.size(chi)):
        result = integrate.quad(lambda x: kv(5./3.,x), y[i], 1e4)
        F_chi[i] = 4*chi[i]**2/eta**2*y[i]*kv(2./3.,y[i])+(1-2*chi[i]/eta)*y[i]*result[0]
#       print(result)
        result = integrate.quad(lambda x: kv(5./3.,x), y_c[i], 1e4)
        F_chi_c[i]=y_c[i]*result[0]
    print(i)
    
    hso_dy[i] = np.sum(F_chi/chi*chi*(10**(0.5*d_log_chi)-10**(-0.5*d_log_chi)))
    hso_dy_c[i] = np.sum(F_chi_c/chi*chi*(10**(0.5*d_log_chi)-10**(-0.5*d_log_chi)))
print('qed:',hso_dy)    
print('classic:',hso_dy_c)    
