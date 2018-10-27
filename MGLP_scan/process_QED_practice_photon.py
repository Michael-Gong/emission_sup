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
part_number=100000
nsteps=0
insert1='./Data_qed/'
insert_n='_0'

px=np.loadtxt(insert1+'px'+insert_n+'.txt')
py=np.loadtxt(insert1+'py'+insert_n+'.txt')
gg=(px**2+py**2+1)**0.5
b0=np.loadtxt(insert1+'bz_part'+insert_n+'.txt')
t0=np.loadtxt(insert1+'t'+insert_n+'.txt')
gg = gg[0]
#b0 = 100 #b0[3]
t0 = t0[-1]
alpha_0 = 0.1
#omega_critic = 1.5*gg**3/r0
omega_critic = 3*gg**2.5*alpha_0**0.5
r0=gg/(2*alpha_0**0.5*gg**0.5)

#circle_number=t0/(2*np.pi*gg/b0)

photon=np.loadtxt(insert1+'qed_photon'+insert_n+'.txt')
photon_py=np.loadtxt(insert1+'qed_py_photon'+insert_n+'.txt')
photon_px=np.loadtxt(insert1+'qed_px_photon'+insert_n+'.txt')
weight_photon=photon
photon = photon*m0*v0**2/(h_planck/6.28*frequency)

theta_grid = np.linspace(-100.5,100.5,202)
theta_bin  = np.linspace(-100,100,201)

def pxpy_to_energy(gamma, weight):
      binsize = 1000
      dd_grid = 6.0/1000.0
      en_grid = np.logspace(np.log10(1e-5*omega_critic)+dd_grid*0.5,np.log10(10*omega_critic)-dd_grid*0.5,binsize)
      en_bin  = np.logspace(np.log10(1e-5*omega_critic),np.log10(10*omega_critic),binsize+1)
      en_value = np.zeros_like(en_grid) 
      for i in range(binsize):
        en_value[i] = np.sum(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1]) ])/(en_bin[i+1]-en_bin[i])
        print('number:',np.size(weight[ (en_bin[i]<=gamma) & (gamma<en_bin[i+1])]),'; value:',en_value[i])
      return (en_grid, en_value)

for i_theta in range(len(theta_bin)):
    theta_photon = np.arctan2(photon_py,photon_px)/np.pi*180
    
    value_x,value_y = pxpy_to_energy(photon[ (theta_grid[i_theta]<theta_photon) & (theta_photon<=theta_grid[i_theta+1])],weight_photon[(theta_grid[i_theta]<theta_photon) & (theta_photon<=theta_grid[i_theta+1])])
    print(value_y)
    
    #print('total circle: ',circle_number)
    
    value_y = value_y/part_number/((theta_grid[-1]-theta_grid[-2])/180*np.pi)*2*np.pi
    
    #theory_line = norm_fac*3*(2*grid_omega_x*gg/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.6667,xi)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.33333,xi)**2)
    
    ratio_c_to_q = np.zeros_like(value_y)
    for i in range(len(ratio_c_to_q)):
        result = (integrate.quad(lambda x: kv(1.6666666667,x), value_x[i]/omega_critic, 1e4))
        ratio_c_to_q[i] = (3)**0.5/4/np.pi**2*(value_x[i]/omega_critic)*gg*(kv(0.66666666666667,value_x[i]/omega_critic/2))**2/result[0]
    
    #norm_x = matplotlib.colors.Normalize()
    plt.subplot(1,2,1)
    #plt.plot(grid_omega_x,theory_line,'-r',linewidth=3,label='theoretical equation')
    plt.plot(value_x, value_y*ratio_c_to_q, markeredgecolor='b',marker='o',linestyle='None',linewidth=3,label='revised QED_photon')
    plt.plot(value_x, value_y, markeredgecolor='k',linewidth=3,marker='o',linestyle='None',label='QED_photon')
    
    norm_fac = q0**2/(4*pi*epsilon0*4*pi**2*v0)/(m0*v0**2)*frequency
    theta_1=0
    xi =  value_x*r0/3.0/gg**3*(1.0+gg**2*theta_1**2)**1.5
    theory_line = norm_fac*3*(2*value_x*r0/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi)**2)
    
    #norm_fac_2 = q0**2*3**0.5/(4*pi*epsilon0*v0)/(m0*v0**2)*frequency
    if (theta_bin[i_theta] >0) and (theta_bin[i_theta] <90): 
         b1 = 2*alpha_0**0.5*((gg**2-1)/((np.tan(theta_bin[i_theta]/pi2d))**2+1))**0.25 
         omega_critic_1 = 1.5*gg**2*b1
         r1 = gg/b1  
         xi_1 =  value_x*r1/3.0/gg**3*(1.0+gg**2*theta_1**2)**1.5
         theory_line_1 = norm_fac*3*(2*value_x*r1/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi_1)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi_1)**2)
    
    plt.plot(value_x,theory_line,':r',linewidth=2,label='theoretical equation') 
    if (theta_bin[i_theta] >0) and (theta_bin[i_theta] <90): 
        plt.plot(value_x,theory_line_1,':b',linewidth=2,label='theoretical equation with r') 
    
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
    plt.xlim(0,10*omega_critic)
    #plt.ylim(0,8.5e-7)
    plt.legend(loc='best',fontsize=16,framealpha=1.0)
    #plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)
    
    plt.subplot(1,2,2)
    #plt.plot(grid_omega_x, data_omega, '-k',linewidth=1)
    #plt.plot(grid_omega_x,theory_line,'-r',linewidth=3)
    #plt.plot(x_line, y_line, ':b',linewidth=3)
    #cbar=plt.colorbar(ticks=np.linspace(0.0, 4, 5))
    #cbar.set_label(r'$log_{10}\frac{dI}{\sin\theta d\theta d\omega}$'+' [A.U.]', fontdict=font)
    #cbar.ax.set_yticklabels(cbar.ax.get_yticklabels(), fontsize=20)
    #plt.plot(x_omega,np.sum(np.sum(data_I_t,axis=0),axis=0),'-b',linewidth=3)
    #### manifesting colorbar, changing label and axis properties ####
    #plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
    plt.plot(value_x, value_y*ratio_c_to_q, markeredgecolor='b',marker='o',linestyle='None',linewidth=3,label='revised QED_photon')
    plt.plot(value_x, value_y, markeredgecolor='k',linewidth=3,marker='o',linestyle='None',label='QED_photon')
    plt.plot(value_x,theory_line,':r',linewidth=2,label='theoretical equation') 
    if (theta_bin[i_theta] >0) and (theta_bin[i_theta] <90): 
        plt.plot(value_x,theory_line_1,':b',linewidth=2,label='theoretical equation with r') 
    plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
    plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
    plt.xticks(fontsize=20); plt.yticks(fontsize=20);
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0,10*omega_critic)
    plt.ylim(1e-8,3e-5)
    #plt.ylim(1e-9,10e-7)
    #plt.legend(loc='upper right',fontsize=16,framealpha=1.0)
    #plt.text(285,4e8,'t='+str(round(time/1.0e-15,0))+' fs',fontdict=font)
    #plt.subplots_adjust(left=0.2, bottom=None, right=0.88, top=None,wspace=None, hspace=None)
    
    fig = plt.gcf()
    fig.set_size_inches(18.0, 8.5)
    fig.savefig(insert1+'spectral_qed_omega_phi='+str(int(theta_bin[i_theta])).zfill(4)+'.png',format='png',dpi=160)
    plt.close("all")

