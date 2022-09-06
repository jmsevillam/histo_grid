import numpy as np
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm

import matplotlib as mpl
mpl.rcParams['mathtext.fontset'] = 'cm'

x,y,z,d=np.genfromtxt('data_hist.dat',skip_header=1).T

fig=plt.figure(figsize=(10,10))
ax=fig.add_subplot(111, projection='3d')
pmf = ax.scatter(x, y, z,c=d,s=10,marker='o')

plt.show() 

