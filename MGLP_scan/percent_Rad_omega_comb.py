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
from decimal import Decimal

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

###############read data into array################
from_path = './C-Rad/'
grid_omega_x  = np.loadtxt(from_path+'grid_omega_x.txt')
grid_theta_y  = np.loadtxt(from_path+'grid_theta_y.txt')
grid_phi_z    = np.loadtxt(from_path+'grid_phi_z.txt')
data_I        = np.loadtxt(from_path+'data.txt')

#print(grid_omega_x.size)
#print(grid_theta_y*pi2d)
#print(grid_phi_z*pi2d)
#print(data_I.size)

data_I  =  data_I.reshape(grid_omega_x.size, grid_theta_y.size, grid_phi_z.size)

#print(data_I.shape)

omega_bin_size = np.zeros_like(grid_omega_x)
lg10d  = np.log10(grid_omega_x[1])-np.log10(grid_omega_x[0])
for i in range(len(omega_bin_size)):
    omega_bin_size[i] = grid_omega_x[i]*(10**(lg10d/2)-10**(-lg10d/2))

px=np.loadtxt('./Data/px_0.txt')
py=np.loadtxt('./Data/py_0.txt')
gg=(px**2+py**2+1)**0.5
b0=np.loadtxt('./Data/bz_part_0.txt')
t0=np.loadtxt('./Data/t_0.txt')
gg = gg[0]
#b0 = 100 #b0[3]
t0 = t0[-1]
#r0=gg/b0
alpha_0 = 1
#omega_critic = 1.5*gg**3/r0
omega_critic = 3*gg**2.5*alpha_0**0.5
r0=gg/(2*alpha_0**0.5*gg**0.5)


for i_theta in range(grid_theta_y.size):
  for i_phi in range(grid_phi_z.size):
    if grid_theta_y.size >1:
        theta_1 = abs(grid_theta_y[i_theta] - pi/2)
    else:
        theta_1 = abs(grid_theta_y - pi/2)
    
#    data_omega = np.sum(np.sum(data_I,axis=2),axis=1)
    data_omega = data_I
    norm_fac = q0**2/(4*pi*epsilon0*4*pi**2*v0)/(m0*v0**2)*frequency
    data_omega = data_omega*norm_fac
    
    xi =  grid_omega_x*r0/3.0/gg**3*(1.0+gg**2*theta_1**2)**1.5
    
    y_line = np.linspace(1e-7,5e-4,1000)
    x_line = np.zeros_like(y_line)+omega_critic
   
    if (abs(grid_phi_z[i_phi]*pi2d) <89.5): 
         b1 = 2*alpha_0**0.5*((gg**2-1)/((np.tan(grid_phi_z[i_phi]))**2+1))**0.25 
         omega_critic_1 = 1.5*gg**2*b1
         r1 = gg/b1  
         xi_1 =  grid_omega_x*r1/3.0/gg**3*(1.0+gg**2*theta_1**2)**1.5
         theory_line_1 = norm_fac*3*(2*grid_omega_x*r1/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi_1)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi_1)**2)
    theory_line = norm_fac*3*(2*grid_omega_x*r0/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi)**2)
    
    if (abs(grid_phi_z[i_phi]*pi2d) <89.5): 
         theory_sum = np.sum((theory_line_1*omega_bin_size)[grid_omega_x < 1000])
         calculate_sum = np.sum((data_omega[:,i_theta,i_phi]*omega_bin_size)[grid_omega_x < 1000])
         #print(grid_theta_y,grid_phi_z)
         print('phi='+str(int(round(grid_phi_z[i_phi]*pi2d,1))).zfill(4)+'; theory:','%.4E' % Decimal(theory_sum),';;; calculate:','%.4E' % Decimal(calculate_sum),';;; ratio:','%.4E' % Decimal(calculate_sum/theory_sum*100),'%')
