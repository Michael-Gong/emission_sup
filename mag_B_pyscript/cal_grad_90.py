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
from_path = './in_90_grad/'
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
b0=2*(alpha*py0*(1-np.cos(90.0/180.0*np.pi)))**0.5*1 #0.1 is for ratio
omega_critic = 1.5*gg**2*b0
print('critical_freq:',omega_critic)
beta=(1.-1./gg**2)**0.5
r0=gg*beta/b0
theta_1 = abs(grid_theta_y - pi/2)
#    data_omega = np.sum(np.sum(data_I,axis=2),axis=1)
data_omega = data_I
norm_fac = q0**2/(4*pi*epsilon0*4*pi**2*v0)/(m0*v0**2)*frequency
data_omega = data_omega*norm_fac
xi =  grid_omega_x*r0/3.0/gg**3*(1.0+gg**2*theta_1**2)**1.5
theory_line = norm_fac*3*(2*grid_omega_x*r0/3.0/gg**2)**2*(1+gg**2*theta_1**2)**2*(kv(0.66666666666667,xi)**2+(gg**2*theta_1**2)/(1.0+gg**2*theta_1**2)*kv(0.3333333333333,xi)**2)


segmengt_list=['1e-4~~3e-4','1e-3~~3e-3','1e-2~~3e-2','1e-1~~3e-1','1e0~~3e0','1e1~~3e1']
seg_list_num =np.array([1e-4,1e-3,1e-2,1e-1,1e0,1e1])
sim_seg_data =np.zeros([7,6])
the_seg_data =np.zeros(6)

theta_list=['theta=90.000','theta=90.090','theta=90.900','theta=99.000','theta=81.000','theta=89.100','theta=89.910']
theta_j   =[0,2,3,4,8,7,6]
delta_grid=grid_omega_x[1:]-grid_omega_x[:-1]

seg_left = np.zeros(6, dtype=np.int16)
seg_right= np.zeros(6, dtype=np.int16)

for i in range(6):
    seg_left[i]  = int(np.min(np.where(grid_omega_x>seg_list_num[i]*omega_critic)))
    seg_right[i] = int(np.min(np.where(grid_omega_x>3*seg_list_num[i]*omega_critic)))
    the_seg_data[i] = np.sum(theory_line[seg_left[i]:seg_right[i]]*(delta_grid[seg_left[i]:seg_right[i]]))

print(seg_left)
print(seg_right)
print(the_seg_data)

for j in range(7):
    for i in range(6):
        sim_seg_data[j,i] = np.sum(data_omega[seg_left[i]:seg_right[i],0,theta_j[j]]*(delta_grid[seg_left[i]:seg_right[i]])) 
    print('theta='+theta_list[j]+':',)
    print(sim_seg_data[j,:])
    print(sim_seg_data[j,:]/the_seg_data)

the_table = plt.table(cellText=sim_seg_data/the_seg_data, rowLabels=theta_list, colLabels=segmengt_list, loc='center')
plt.axis('off')

fig = plt.gcf()
fig.set_size_inches(20.0, 16.0)
fig.savefig('./table_grad_90.png',format='png',dpi=160)
plt.close("all")
    
#plt.subplot(1,3,1)
#plt.plot(grid_omega_x,theory_line,'--k',linewidth=3,label='Synchrotron radiation')
#plt.plot(grid_omega_x, data_omega[:,0,0],linestyle='-',color='crimson',linewidth=3,label=r'$\theta_d=15.000^\circ$')
#plt.plot(grid_omega_x, data_omega[:,0,2],linestyle='-',color='darkorange',linewidth=3,label=r'$\theta_d=15.015^\circ$')
#plt.plot(grid_omega_x, data_omega[:,0,3],linestyle='-',color='mediumseagreen',linewidth=3,label=r'$\theta_d=15.150^\circ$')
#plt.plot(grid_omega_x, data_omega[:,0,4],linestyle='-',color='royalblue',linewidth=3,label=r'$\theta_d=16.500^\circ$')
#
#plt.subplot(1,3,2)
#plt.plot(grid_omega_x,theory_line,'--k',linewidth=3,label='Synchrotron radiation')
#plt.plot(grid_omega_x, data_omega[:,0,8],linestyle='-',color='royalblue',linewidth=3,label=r'$\theta_d=13.500^\circ$')
#plt.plot(grid_omega_x, data_omega[:,0,7],linestyle='-',color='mediumseagreen',linewidth=3,label=r'$\theta_d=14.850^\circ$')
#plt.plot(grid_omega_x, data_omega[:,0,6],linestyle='-',color='darkorange',linewidth=3,label=r'$\theta_d=14.985^\circ$')
#plt.plot(grid_omega_x, data_omega[:,0,0],linestyle='-',color='crimson',linewidth=3,label=r'$\theta_d=15.000^\circ$')




