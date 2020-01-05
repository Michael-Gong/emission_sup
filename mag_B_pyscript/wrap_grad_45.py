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


font_size =28
###############read data into array################
from_path = './in_45_grad/'
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

py0=100.0
gg= (py0**2+1)**0.5
alpha=10.0
b0=2*(alpha*py0*(1-np.cos(45.0/180.0*np.pi)))**0.5*1 #0.1 is for ratio
#b0=10. #np.loadtxt(from_path+'bz_part_0000.txt')

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
 
plt.subplot(1,3,1)
#plt.plot(grid_omega_x, data_omega[:,0,6], '-k',marker='o',linewidth=2)
plt.plot(grid_omega_x,theory_line,'--k',linewidth=3,label='Synchrotron radiation')
plt.plot(grid_omega_x, data_omega[:,0,0],linestyle='-',color='crimson',linewidth=3,label=r'$\theta_d=45.000^\circ$')
plt.plot(grid_omega_x, data_omega[:,0,2],linestyle='-',color='darkorange',linewidth=3,label=r'$\theta_d=45.045^\circ$')
plt.plot(grid_omega_x, data_omega[:,0,3],linestyle='-',color='mediumseagreen',linewidth=3,label=r'$\theta_d=45.450^\circ$')
plt.plot(grid_omega_x, data_omega[:,0,4],linestyle='-',color='royalblue',linewidth=3,label=r'$\theta_d=49.500^\circ$')
#plt.plot(xx, yy0*norm_fac2, linestyle='-',color='crimson',linewidth=3,label=r'$\theta_d=90.00^\circ$')
#plt.plot(xx, yy1*norm_fac2, linestyle='-',color='darkorange',linewidth=3,label=r'$\theta_d=90.09^\circ$')
#plt.plot(xx, yy2*norm_fac2, linestyle='-',color='mediumseagreen',linewidth=3,label=r'$\theta_d=90.9^\circ$')
#plt.plot(xx, yy3*norm_fac2, linestyle='-',color='royalblue',linewidth=3,label=r'$\theta_d=99.0^\circ$')
#plt.plot(xx, yy4*norm_fac2, linestyle='-',color='purple',linewidth=2,label=r'$\theta_d=108.0^\circ$')
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
plt.xticks(fontsize=font_size); 
plt.yticks(fontsize=font_size);
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5*omega_critic,5e1*omega_critic)
plt.ylim(1e-9,1e-4)
plt.legend(loc='best',fontsize=font_size-8,framealpha=0.5)
plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)


plt.subplot(1,3,2)
#plt.plot(grid_omega_x, data_omega[:,0,6], '-k',marker='o',linewidth=2)
plt.plot(grid_omega_x,theory_line,'--k',linewidth=3,label='Synchrotron radiation')
plt.plot(grid_omega_x, data_omega[:,0,8],linestyle='-',color='royalblue',linewidth=3,label=r'$\theta_d=40.500^\circ$')
plt.plot(grid_omega_x, data_omega[:,0,7],linestyle='-',color='mediumseagreen',linewidth=3,label=r'$\theta_d=44.550^\circ$')
plt.plot(grid_omega_x, data_omega[:,0,6],linestyle='-',color='darkorange',linewidth=3,label=r'$\theta_d=44.955^\circ$')
plt.plot(grid_omega_x, data_omega[:,0,0],linestyle='-',color='crimson',linewidth=3,label=r'$\theta_d=45.000^\circ$')
#plt.plot(xx, yy0*norm_fac2, linestyle='-',color='crimson',linewidth=3,label=r'$\theta_d=90.00^\circ$')
#plt.plot(xx, yy5*norm_fac2, linestyle='-',color='darkorange',linewidth=3,label=r'$\theta_d=89.91^\circ$')
#plt.plot(xx, yy6*norm_fac2, linestyle='-',color='mediumseagreen',linewidth=3,label=r'$\theta_d=89.1^\circ$')
#plt.plot(xx, yy7*norm_fac2, linestyle=':',color='skyblue',linewidth=2)
#plt.plot(xx, yy8*norm_fac2, linestyle=':',color='purple',linewidth=2)
plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
#plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
plt.xticks(fontsize=font_size); 
plt.yticks(fontsize=font_size);
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5*omega_critic,5e1*omega_critic)
plt.ylim(1e-9,1e-4)
plt.legend(loc='best',fontsize=font_size-8,framealpha=0.5)
plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)



plt.subplot(1,3,3)
insert1='./qe_45_grad/'
px=70.71 # np.loadtxt(insert1+'px_0000.txt')
py=70.71 #np.loadtxt(insert1+'py_0000.txt')
p0=(px**2+py**2)**0.5
gg=(px**2+py**2+1)**0.5
alpha=10.0
b0=2*(alpha*p0*(1-np.cos(45.0/180.0*np.pi)))**0.5*1 #0.1 is for ratio
insert2='./in_45_grad/'
to_px=np.loadtxt(insert2+'px_0000.txt')
to_py=np.loadtxt(insert2+'py_0000.txt')
to_theta = np.arctan2(to_py,to_px)/np.pi*180.

t0 = 14.611926
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
grid_x = ppx*4*omega_critic/3/eta**2
grid_y = ppy*3*eta**2/4/omega_critic
xi =  grid_x*r0/3.0/gg**3
for i in range(np.size(grid_x)):
    result = integrate.quad(lambda x: kv(1.6666666667,x), 2*xi[i], 1e3*2*xi[i])
    grid_y[i] = grid_y[i]*(2*xi[i]*gg*3**0.5/(4*pi**2)*kv(0.66666666666667,xi[i])**2)/result[0]
plt.bar(grid_x, grid_y, width=0.3*grid_x, color='blue',alpha=0.5,edgecolor='black',linewidth=2,label='Monte-carlo method')
#plt.plot(grid_x,grid_y,marker='o',markerfacecolor='tomato',markeredgecolor='darkred',markersize=6.6,markeredgewidth=2,label=r'$\theta_d=15.000^\circ$',linestyle='None')

to_t=np.loadtxt(insert2+'t_0000.txt')
to_t=to_t[0:np.size(to_t)//2]
to_theta=to_theta[0:np.size(to_theta)//2]

print(to_t)
to_t_2=to_t[(to_theta>-46)&(to_theta<-44)]
print(to_t_2)
t0=(to_t_2[-1]-to_t_2[0])*2.0 #14.11078 #np.loadtxt(insert1+'t_0000.txt')
insert_n='_0000'
photon_px=np.loadtxt(insert1+'photon_px'+insert_n+'.txt')
photon_py=np.loadtxt(insert1+'photon_py'+insert_n+'.txt')
photon_theta=np.arctan2(photon_py,photon_px)/np.pi*180.
photon_px=photon_px[(photon_theta>-46)&(photon_theta<-44)]
photon_py=photon_py[(photon_theta>-46)&(photon_theta<-44)]
photon = (photon_px**2+photon_py**2)**0.5
chi   = photon*b0/2.0/4.1e5
eta   = gg*b0/4.1e5
dt    = t0/6.28
r0=gg/b0
omega_critic=1.5*gg**3/r0
weight_photon = photon/dt/1e4
ppx,ppy = pxpy_to_energy_rr(chi, weight_photon, eta)
grid_x = ppx*4*omega_critic/3/eta**2
grid_y = ppy*3*eta**2/4/omega_critic
xi =  grid_x*r0/3.0/gg**3
for i in range(np.size(grid_x)):
    result = integrate.quad(lambda x: kv(1.6666666667,x), 2*xi[i], 1e3*2*xi[i])
    grid_y[i] = grid_y[i]*(2*xi[i]*gg*3**0.5/(4*pi**2)*kv(0.66666666666667,xi[i])**2)/result[0]
#plt.bar(grid_x, grid_y, width=0.3*grid_x, color='blue',alpha=0.5,edgecolor='black',linewidth=2,label='Monte-carlo method')
plt.plot(grid_x,grid_y,marker='o',markerfacecolor='tomato',markeredgecolor='darkred',markersize=10,markeredgewidth=3,label=r'$\theta_d=45.000^\circ$',linestyle='None')

print(to_t)
to_t_2=to_t[(to_theta>-41.5)&(to_theta<-39.5)]
print(to_t_2)
t0=(to_t_2[-1]-to_t_2[0])*2.0 #14.11078 #np.loadtxt(insert1+'t_0000.txt')
insert_n='_0000'
photon_px=np.loadtxt(insert1+'photon_px'+insert_n+'.txt')
photon_py=np.loadtxt(insert1+'photon_py'+insert_n+'.txt')
photon_theta=np.arctan2(photon_py,photon_px)/np.pi*180.
photon_px=photon_px[(photon_theta>-41.5)&(photon_theta<-39.5)]
photon_py=photon_py[(photon_theta>-41.5)&(photon_theta<-39.5)]
photon = (photon_px**2+photon_py**2)**0.5
chi   = photon*b0/2.0/4.1e5
eta   = gg*b0/4.1e5
dt    = t0/6.28
r0=gg/b0
omega_critic=1.5*gg**3/r0
weight_photon = photon/dt/1e4
ppx,ppy = pxpy_to_energy_rr(chi, weight_photon, eta)
grid_x = ppx*4*omega_critic/3/eta**2
grid_y = ppy*3*eta**2/4/omega_critic
xi =  grid_x*r0/3.0/gg**3
for i in range(np.size(grid_x)):
    result = integrate.quad(lambda x: kv(1.6666666667,x), 2*xi[i], 1e3*2*xi[i])
    grid_y[i] = grid_y[i]*(2*xi[i]*gg*3**0.5/(4*pi**2)*kv(0.66666666666667,xi[i])**2)/result[0]
#plt.bar(grid_x, grid_y, width=0.3*grid_x, color='blue',alpha=0.5,edgecolor='black',linewidth=2,label='Monte-carlo method')
plt.plot(grid_x,grid_y,marker='o',markerfacecolor='skyblue',markeredgecolor='darkblue',markersize=10,markeredgewidth=3,label=r'$\theta_d=40.500^\circ$',linestyle='None')

print(to_t)
to_t_2=to_t[(to_theta>-31)&(to_theta<-29)]
print(to_t_2)
t0=(to_t_2[-1]-to_t_2[0])*2.0 #14.11078 #np.loadtxt(insert1+'t_0000.txt')
insert_n='_0000'
photon_px=np.loadtxt(insert1+'photon_px'+insert_n+'.txt')
photon_py=np.loadtxt(insert1+'photon_py'+insert_n+'.txt')
photon_theta=np.arctan2(photon_py,photon_px)/np.pi*180.
photon_px=photon_px[(photon_theta>-31)&(photon_theta<-29)]
photon_py=photon_py[(photon_theta>-31)&(photon_theta<-29)]
photon = (photon_px**2+photon_py**2)**0.5
chi   = photon*b0/2.0/4.1e5
eta   = gg*b0/4.1e5
dt    = t0/6.28
r0=gg/b0
omega_critic=1.5*gg**3/r0
weight_photon = photon/dt/1e4
ppx,ppy = pxpy_to_energy_rr(chi, weight_photon, eta)
grid_x = ppx*4*omega_critic/3/eta**2
grid_y = ppy*3*eta**2/4/omega_critic
xi =  grid_x*r0/3.0/gg**3
for i in range(np.size(grid_x)):
    result = integrate.quad(lambda x: kv(1.6666666667,x), 2*xi[i], 1e3*2*xi[i])
    grid_y[i] = grid_y[i]*(2*xi[i]*gg*3**0.5/(4*pi**2)*kv(0.66666666666667,xi[i])**2)/result[0]
#plt.bar(grid_x, grid_y, width=0.3*grid_x, color='blue',alpha=0.5,edgecolor='black',linewidth=2,label='Monte-carlo method')
plt.plot(grid_x,grid_y,marker='o',markerfacecolor='yellow',markeredgecolor='darkorange',markersize=10,markeredgewidth=3,label=r'$\theta_d=30.000^\circ$',linestyle='None')

print(to_t)
to_t_2=to_t[(to_theta>-1)&(to_theta<1)]
print(to_t_2)
t0=to_t_2[-1]-to_t_2[0] #14.11078 #np.loadtxt(insert1+'t_0000.txt')
insert_n='_0000'
photon_px=np.loadtxt(insert1+'photon_px'+insert_n+'.txt')
photon_py=np.loadtxt(insert1+'photon_py'+insert_n+'.txt')
photon_theta=np.arctan2(photon_py,photon_px)/np.pi*180.
photon_px=photon_px[(photon_theta>-1)&(photon_theta<1)]
photon_py=photon_py[(photon_theta>-1)&(photon_theta<1)]
photon = (photon_px**2+photon_py**2)**0.5
chi   = photon*b0/2.0/4.1e5
eta   = gg*b0/4.1e5
dt    = t0/6.28
r0=gg/b0
omega_critic=1.5*gg**3/r0
weight_photon = photon/dt/1e4
ppx,ppy = pxpy_to_energy_rr(chi, weight_photon, eta)
grid_x = ppx*4*omega_critic/3/eta**2
grid_y = ppy*3*eta**2/4/omega_critic
xi =  grid_x*r0/3.0/gg**3
for i in range(np.size(grid_x)):
    result = integrate.quad(lambda x: kv(1.6666666667,x), 2*xi[i], 1e3*2*xi[i])
    grid_y[i] = grid_y[i]*(2*xi[i]*gg*3**0.5/(4*pi**2)*kv(0.66666666666667,xi[i])**2)/result[0]
#plt.bar(grid_x, grid_y, width=0.3*grid_x, color='blue',alpha=0.5,edgecolor='black',linewidth=2,label='Monte-carlo method')
plt.plot(grid_x,grid_y,marker='o',markerfacecolor='magenta',markeredgecolor='black',markersize=10,markeredgewidth=3,label=r'$\theta_d=0.000^\circ$',linestyle='None')



grid_omega_x = np.logspace(0,np.log10(100*omega_critic),1000)
norm_fac = q0**2/(4*pi*epsilon0*4*pi**2*v0)/(m0*v0**2)*frequency
theta_1 = 0
xi =  grid_omega_x*r0/3.0/gg**3*(1.0+gg**2*theta_1**2)**1.5
theory_line = norm_fac*3*(2*grid_omega_x*r0/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi)**2)
theory_line_1 = np.zeros_like(theory_line)
for i in range(np.size(xi)):
    result = integrate.quad(lambda x: kv(1.6666666667,x), 2*xi[i], 1e3*2*xi[i])
    theory_line_1[i] = theory_line[i]/(2*xi[i]*gg*3**0.5/(4*pi**2)*kv(0.66666666666667,xi[i])**2)*result[0]
plt.plot(grid_omega_x,theory_line,':',color='lime',linewidth=6,label='Synchrotron radiation')

plt.xlabel('$\omega$ [$\omega_0$]',fontdict=font)
#plt.ylabel(r'$\frac{dI^2}{d\omega d\Omega}$'+' [$m_ec^2/\omega_0$]',fontdict=font)
plt.xticks(fontsize=font_size); 
plt.yticks(fontsize=font_size);
plt.xscale('log')
plt.yscale('log')
plt.xlim(1e-5*omega_critic,5e1*omega_critic)
plt.ylim(1e-9,1e-4)
plt.legend(loc='best',fontsize=font_size-8,framealpha=0.5)
plt.grid(which='major',color='k', linestyle='--', linewidth=0.3)
#plt.grid(which='minor',color='k', linestyle=':', linewidth=0.1)



plt.subplots_adjust(left=0.15, bottom=0.12, right=0.99, top=0.97, wspace=0.18, hspace=0.2) 
fig = plt.gcf()
fig.set_size_inches(28.5, 9.0)
fig.savefig(to_path+'wrap_grad_in45.png',format='png',dpi=160)
#    fig.savefig(to_path+'spectral_theta='+str(int(np.round(grid_theta_y*pi2d,1))).zfill(4)+'_phi='+str(int(i_phi)).zfill(4)+'.png',format='png',dpi=160)
plt.close("all")
    
