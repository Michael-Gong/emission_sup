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

font_size =28
###############read data into array################
from_path = './Data_full/'
to_path = './'
grid_omega_x  = np.loadtxt(from_path+'grid_omega_x.txt')
grid_theta_y  = np.loadtxt(from_path+'grid_theta_y.txt')
grid_phi_z    = np.loadtxt(from_path+'grid_phi_z.txt')
data_I        = np.loadtxt(from_path+'grid_data.txt')

print(grid_omega_x.size)
print(grid_theta_y*pi2d)
print(grid_phi_z*pi2d)
print(data_I.size)

data_I  =  data_I.reshape(grid_omega_x.size, grid_theta_y.size, grid_phi_z.size)

print(data_I.shape)

gg= (100.**2+1)**0.5
b0=10. #np.loadtxt(from_path+'bz_part_0000.txt')
omega_critic = 1.5*gg**2*b0
beta=(1.-1./gg**2)**0.5
r0=gg*beta/b0

theta_1 = abs(grid_theta_y - pi/2)
    
#    data_omega = np.sum(np.sum(data_I,axis=2),axis=1)
data_omega = data_I
norm_fac = q0**2/(4*pi*epsilon0*4*pi**2*v0)/(m0*v0**2)*frequency
data_omega = data_omega*norm_fac
    
xi =  grid_omega_x*r0/3.0/gg**3*(1.0+gg**2*theta_1**2)**1.5
    
y_line = np.linspace(1e-7,5e-4,1000)
x_line = np.zeros_like(y_line)+omega_critic
    
theory_line = norm_fac*3*(2*grid_omega_x*r0/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi)**2)

from_path2='./analytical_2/'
xx  = np.loadtxt(from_path2+'data_grid.txt')    
yy0 = np.loadtxt(from_path2+'data_0.txt')    
yy1 = np.loadtxt(from_path2+'data_1.txt')    
yy2 = np.loadtxt(from_path2+'data_2.txt')    
yy3 = np.loadtxt(from_path2+'data_3.txt')    
yy4 = np.loadtxt(from_path2+'data_4.txt')    
yy5 = np.loadtxt(from_path2+'data_5.txt')    
yy6 = np.loadtxt(from_path2+'data_6.txt')    
yy7 = np.loadtxt(from_path2+'data_7.txt')    
yy8 = np.loadtxt(from_path2+'data_8.txt')    

norm_fac2=4*beta**2*norm_fac
 
plt.subplot(2,2,2)
#plt.plot(grid_omega_x, data_omega[:,0,6], '-k',marker='o',linewidth=2)
plt.plot(grid_omega_x,theory_line,'--k',linewidth=3,label='Synchrotron radiation')
plt.plot(xx, yy0*norm_fac2, linestyle='-',color='crimson',linewidth=3,label=r'$\theta_d=90.00^\circ$')
plt.plot(xx, yy1*norm_fac2, linestyle='-',color='darkorange',linewidth=3,label=r'$\theta_d=90.09^\circ$')
plt.plot(xx, yy2*norm_fac2, linestyle='-',color='mediumseagreen',linewidth=3,label=r'$\theta_d=90.9^\circ$')
plt.plot(xx, yy3*norm_fac2, linestyle='-',color='royalblue',linewidth=3,label=r'$\theta_d=99.0^\circ$')
#plt.plot(xx, yy4*norm_fac2, linestyle='-',color='purple',linewidth=2,label=r'$\theta_d=108.0^\circ$')
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
plt.xticks(fontsize=font_size); 
plt.yticks(fontsize=font_size);
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5*omega_critic,1e2*omega_critic)
plt.ylim(1e-10,1e-4)
plt.legend(loc='best',fontsize=font_size-5,framealpha=0.5)
plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)


plt.subplot(2,2,3)
#plt.plot(grid_omega_x, data_omega[:,0,6], '-k',marker='o',linewidth=2)
plt.plot(grid_omega_x,theory_line,'--k',linewidth=3,label='Synchrotron radiation')
plt.plot(xx, yy0*norm_fac2, linestyle='-',color='crimson',linewidth=3,label=r'$\theta_d=90.00^\circ$')
plt.plot(xx, yy5*norm_fac2, linestyle='-',color='darkorange',linewidth=3,label=r'$\theta_d=89.91^\circ$')
plt.plot(xx, yy6*norm_fac2, linestyle='-',color='mediumseagreen',linewidth=3,label=r'$\theta_d=89.1^\circ$')
#plt.plot(xx, yy7*norm_fac2, linestyle=':',color='skyblue',linewidth=2)
#plt.plot(xx, yy8*norm_fac2, linestyle=':',color='purple',linewidth=2)
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
plt.xticks(fontsize=font_size); 
plt.yticks(fontsize=font_size);
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5*omega_critic,1e2*omega_critic)
plt.ylim(1e-10,1e-4)
plt.legend(loc='best',fontsize=font_size-5,framealpha=0.5)
plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)


plt.subplot(2,2,4)
#plt.plot(grid_omega_x, data_omega[:,0,6], '-k',marker='o',linewidth=2)
plt.plot(grid_omega_x,theory_line,'--k',linewidth=3,label='Synchrotron radiation')
#plt.plot(xx, yy0*norm_fac2, linestyle='-',color='crimson',linewidth=3,label=r'$\theta_d=90.00^\circ$')
#plt.plot(xx, yy5*norm_fac2, linestyle='-',color='darkorange',linewidth=3,label=r'$\theta_d=89.91^\circ$')
#plt.plot(xx, yy6*norm_fac2, linestyle='-',color='mediumseagreen',linewidth=3,label=r'$\theta_d=89.1^\circ$')
plt.plot(xx, yy7*norm_fac2, linestyle='-',color='royalblue',linewidth=2,label=r'$\theta_d=81.0^\circ$')
#plt.plot(xx, yy8*norm_fac2, linestyle=':',color='purple',linewidth=2)
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
#plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
plt.xticks(fontsize=font_size); 
plt.yticks(fontsize=font_size);
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5*omega_critic,1e2*omega_critic)
plt.ylim(1e-10,1e-4)
plt.legend(loc='best',fontsize=font_size-5,framealpha=0.5)
plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)


plt.subplots_adjust(left=0.15, bottom=0.12, right=0.99, top=0.98, wspace=0.24, hspace=0.2) 
fig = plt.gcf()
fig.set_size_inches(19.0, 18.0)
fig.savefig(to_path+'wrap_full_con_2.png',format='png',dpi=160)
#    fig.savefig(to_path+'spectral_theta='+str(int(np.round(grid_theta_y*pi2d,1))).zfill(4)+'_phi='+str(int(i_phi)).zfill(4)+'.png',format='png',dpi=160)
plt.close("all")
    
