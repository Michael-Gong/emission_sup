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

###############read data into array################
from_path = './Data_zoomin/'
to_path = from_path
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

px=np.loadtxt(from_path+'px_0000.txt')
py=np.loadtxt(from_path+'py_0000.txt')
gg= (px**2+py**2+1)**0.5
b0=np.loadtxt(from_path+'bz_part_0000.txt')
gg = gg[0]
b0 = abs(np.min(b0)) # b0[1]
r0=gg/b0
omega_critic = 1.5*gg**3/r0

for i_theta in range(grid_theta_y.size):
#for i_theta in range(4,6,1):
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
    
    theory_line = norm_fac*3*(2*grid_omega_x*r0/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi)**2)
    
    #norm_x = matplotlib.colors.Normalize()
    plt.subplot(2,2,1)
    plt.plot(grid_omega_x, data_omega[:,i_theta,i_phi], '-k', marker='o',linewidth=1,label='my code calculation')
    plt.plot(grid_omega_x,theory_line,':r',linewidth=2,label='theoretical equation')
    #plt.plot(x_line, y_line, ':b',linewidth=3,label='$\omega_{c}=(3c\gamma^3)/(2r_0)$')
    #cbar=plt.colorbar(ticks=np.linspace(0.0, 4, 5))
    #cbar.set_label(r'$log_{10}\frac{dI}{\sin\theta d\theta d\omega}$'+' [A.U.]', fontdict=font)
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    #plt.plot(x_omega,np.sum(np.sum(data_I_t,axis=0),axis=0),'-b',linewidth=3)
    #### manifesting colorbar, changing label and axis properties ####
    #plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
    plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
    plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.yscale('log')
#    plt.xlim(2,5990)
#    plt.ylim(0,1e-4)
    plt.legend(loc='best',fontsize=16,framealpha=1.0)
    #plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)
    
    plt.subplot(2,2,2)
    plt.plot(grid_omega_x, data_omega[:,i_theta,i_phi], '-k',marker='o',linewidth=1)
    plt.plot(grid_omega_x,theory_line,':r',linewidth=2)
    #plt.plot(x_line, y_line, ':b',linewidth=3)
    #cbar=plt.colorbar(ticks=np.linspace(0.0, 4, 5))
    #cbar.set_label(r'$log_{10}\frac{dI}{\sin\theta d\theta d\omega}$'+' [A.U.]', fontdict=font)
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    #plt.plot(x_omega,np.sum(np.sum(data_I_t,axis=0),axis=0),'-b',linewidth=3)
    #### manifesting colorbar, changing label and axis properties ####
    #plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
    plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
    plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-5*omega_critic,10*omega_critic)
    plt.ylim(1e-10,1e-4)
    #plt.legend(loc='upper right',fontsize=16,framealpha=1.0)
    #plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)
    
    theta_c = 1/gg*(2*omega_critic/grid_omega_x)**0.33333333333333/np.pi*180
    #norm_x = matplotlib.colors.Normalize()
    plt.subplot(2,2,3)
    plt.plot(grid_omega_x, theta_c, '-k', marker='^',linewidth=1,label='my code calculation')
    #plt.plot(x_line, y_line, ':b',linewidth=3,label='$\omega_{c}=(3c\gamma^3)/(2r_0)$')
    #cbar=plt.colorbar(ticks=np.linspace(0.0, 4, 5))
    #cbar.set_label(r'$log_{10}\frac{dI}{\sin\theta d\theta d\omega}$'+' [A.U.]', fontdict=font)
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    #plt.plot(x_omega,np.sum(np.sum(data_I_t,axis=0),axis=0),'-b',linewidth=3)
    #### manifesting colorbar, changing label and axis properties ####
    #plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
    plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
    plt.ylabel(r'$\theta_c\ [^o]$',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    #plt.yscale('log')
#    plt.xlim(2,5990)
#    plt.ylim(0,1e-4)
    plt.legend(loc='best',fontsize=16,framealpha=1.0)
    plt.xscale('log')
    #plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)
    
    plt.subplot(2,2,4)
    plt.plot(grid_omega_x*1.24e-6, data_omega[:,i_theta,i_phi]*0.51e-13/1.24e-6/(grid_omega_x*1.24e-6*1.6e-19), '-k',marker='o',linewidth=1)
    plt.plot(grid_omega_x*1.24e-6,theory_line*0.51e-13/1.24e-6/(grid_omega_x*1.24e-6*1.6e-19),':r',linewidth=2)
    plt.xlabel(r'$\varepsilon_\gamma$'+' [MeV]',fontdict=font)
    plt.ylabel(r'$\frac{dN^2}{d\varepsilon_\gamma d\Omega}$'+' [MeV$^{-1}$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1e-5*omega_critic*1.24e-6,10*omega_critic*1.24e-6)
    plt.ylim(1e3,1e9)
    #plt.legend(loc='upper right',fontsize=16,framealpha=1.0)
    #plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)
 #   plt.xlim(2,6000)
 #   plt.ylim(1e-7,5e-4)
    #plt.legend(loc='upper right',fontsize=16,framealpha=1.0)
    #plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)
    
    
    
    fig = plt.gcf()
    fig.set_size_inches(18.0, 18.0)
    fig.savefig(to_path+'spectral_theta='+str(int(np.round(grid_theta_y*pi2d,1))).zfill(4)+'_phi='+str(int(round(grid_phi_z[i_phi]*pi2d,1))).zfill(4)+'.png',format='png',dpi=160)
#    fig.savefig(to_path+'spectral_theta='+str(int(np.round(grid_theta_y*pi2d,1))).zfill(4)+'_phi='+str(int(i_phi)).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")
    
