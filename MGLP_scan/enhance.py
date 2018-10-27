#%matplotlib inline
#import sdf
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma
import matplotlib as mpl
mpl.style.use('https://raw.githubusercontent.com/Michael-Gong/DLA_project/master/style')
from matplotlib import colors, ticker, cm
from matplotlib.mlab import bivariate_normal
from optparse import OptionParser
import os
import os
from mpl_toolkits.mplot3d import Axes3D
import random
from mpl_toolkits import mplot3d
from matplotlib import rc
import matplotlib.transforms as mtransforms
#from colour import Color

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

font = {'family' : 'helvetica',
        'color'  : 'black',
		'weight' : 'normal',
        'size'   : 25,
	   }

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

from matplotlib.colors import Normalize
class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

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
denunit    =     frequency**2*epsilon0*m0/q0**2
#print 'electric field unit: '+str(exunit)
#print 'magnetic field unit: '+str(bxunit)
#print 'density unit nc: '+str(denunit)

from_path='./result/'

enhancement=np.load(from_path+'enhance.npy')
axis_b     =np.load(from_path+'axis_b.npy')
axis_p     =np.load(from_path+'axis_p.npy')
axis_phi   =np.load(from_path+'axis_phi.npy')
axis_omega =np.load(from_path+'axis_omega.npy')

for i_phi in range(axis_phi.size):
    for i_omega in range(axis_omega.size):
        data = enhancement[:,:,i_phi,i_omega]
        
        levels = np.logspace(0, 0.4, 40)
        #levels = np.linspace(1, 10, 40)
#plt.contourf(x, y, ex, levels=levels, cmap=cm.pink_r, antialiased=False)
#plt.pcolormesh(x, y, ex, norm=colors.SymLogNorm(linthresh=0.03, linscale=0.03, vmin=np.min(ex.T), vmax=np.max(ex.T)), cmap=cm.pink_r)
        norm = MidpointNormalize(midpoint=50)
        plt.pcolormesh(axis_b, axis_p, data.T, norm=mpl.colors.LogNorm(vmin=min(levels),vmax=max(levels)), cmap=cm.cubehelix_r)
        #plt.pcolormesh(axis_b, axis_p, data.T, norm=mpl.colors.Normalize(vmin=min(levels),vmax=max(levels)), cmap=cm.pink_r)

        plt.axis([axis_b.min(), axis_b.max(), axis_p.min(), axis_p.max()])
        #### manifesting colorbar, changing label and axis properties ####
        cbar=plt.colorbar()#ticks=[np.min(ex), -eee/2, 0, eee/2, np.min()])
        cbar.set_label(r'$\frac{calculation}{theory}$',fontdict=font)        

        plt.xlabel(r'$lg\alpha_0$',fontdict=font)
        plt.ylabel(r'$P_{0}\ [m_ec]$',fontdict=font)
        plt.xticks(fontsize=25); plt.yticks(fontsize=25);
        #plt.title('calculation/theory for '+'phi='+str(int(axis_phi[i_phi]))+'--omega<'+str(int(axis_omega[i_omega])))
        plt.xlim(axis_b.min(), axis_b.max())
        plt.ylim(axis_p.min(), axis_p.max())

        fig = plt.gcf()
        fig.set_size_inches(8, 6.5)
        #fig.set_size_inches(5, 4.5)
        fig.savefig(from_path+'phi='+str(int(axis_phi[i_phi]))+'_omega<'+str(int(axis_omega[i_omega]))+'.png',format='png',dpi=160)
        plt.close("all")
