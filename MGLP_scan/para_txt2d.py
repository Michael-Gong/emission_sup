#%matplotlib inline
#import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import scipy.integrate as integrate
import scipy.special as special 
from scipy.special import kv
from decimal import Decimal


#from colour import Color

######## Constant defined here ########
pi        =     3.1415926535897932384626
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
denunit   =     frequency**2*epsilon0*m0/q0**2
#print 'electric field unit: '+str(exunit)
#print 'magnetic field unit: '+str(bxunit)
#print 'density unit nc: '+str(denunit)

font = {'family' : 'monospace',  
        'color'  : 'black',  
	'weight' : 'normal',  
        'size'   : 20,  
       }  



axis_b=np.linspace(0,40000,41)
#axis_a=np.array([2000,5000,10000,15000])
#axis_a=np.linspace(100,2500,41)
axis_p=np.linspace(10,250,41)
#axis_w=np.linspace(30,90,16)
axis_phi = np.linspace(80,89,10)
axis_omega = np.logspace(3,7,5)
enhancement=np.zeros([41,41,10,6])
#rebo2=np.zeros([16,16])
insert1='./'
insert2='/C-Rad/'
#insert2='epoch2dqe/'
insert_n='_0'
for ib in range(axis_b.size):
    for ip in range(axis_p.size):
        grid_omega_x  = np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+insert2+'grid_omega_x.txt')
        grid_theta_y  = np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+insert2+'grid_theta_y.txt')
        grid_phi_z    = np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+insert2+'grid_phi_z.txt')
        data_I        = np.loadtxt(insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+insert2+'data.txt')

        data_I  =  data_I.reshape(grid_omega_x.size, grid_theta_y.size, grid_phi_z.size)

        #print(data_I.shape)

        omega_bin_size = np.zeros_like(grid_omega_x)
        lg10d  = np.log10(grid_omega_x[1])-np.log10(grid_omega_x[0])
        for i in range(len(omega_bin_size)):
            omega_bin_size[i] = grid_omega_x[i]*(10**(lg10d/2)-10**(-lg10d/2))

        py = axis_p[ip]; px = 0.0
        gg = (px**2+py**2+1)**0.5
        alpha_0 = 10**(axis_b[ib]/1.e4-2.5)

        for i_theta in range(grid_theta_y.size):
            for i_phi in range(grid_phi_z.size):
                if grid_theta_y.size >1:
                    theta_1 = abs(grid_theta_y[i_theta] - pi/2)
                else:
                    theta_1 = abs(grid_theta_y - pi/2)
    
#               data_omega = np.sum(np.sum(data_I,axis=2),axis=1)
                data_omega = data_I
                norm_fac = q0**2/(4*pi*epsilon0*4*pi**2*v0)/(m0*v0**2)*frequency
                data_omega = data_omega*norm_fac
                b1 = 2*alpha_0**0.5*((gg**2-1)/((np.tan(grid_phi_z[i_phi]))**2+1))**0.25 
                omega_critic_1 = 1.5*gg**2*b1
                r1 = gg/b1

                xi_1 =  grid_omega_x*r1/3.0/gg**3*(1.0+gg**2*theta_1**2)**1.5
                theory_line_1 = norm_fac*3*(2*grid_omega_x*r1/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi_1)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi_1)**2)
                for i_omega in range(axis_omega.size):
                    theory_sum = np.sum((theory_line_1*omega_bin_size)[grid_omega_x < axis_omega[i_omega]])
                    calculate_sum = np.sum((data_omega[:,i_theta,i_phi]*omega_bin_size)[grid_omega_x < axis_omega[i_omega]])
                    enhancement[ib,ip,i_phi,i_omega] = calculate_sum/theory_sum

        print("finished "+insert1+'Datab'+str(int(axis_b[ib]))+'p'+str(int(axis_p[ip]))+insert2)

to_path='./result/'

np.save(to_path+'enhance', enhancement)
np.save(to_path+'axis_b', axis_b/1.e4-2.5)
np.save(to_path+'axis_p', axis_p)
np.save(to_path+'axis_phi', axis_phi)
np.save(to_path+'axis_omega', axis_omega)
