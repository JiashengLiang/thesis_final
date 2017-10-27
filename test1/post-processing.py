import numpy as np
import matplotlib.pyplot as plt



"""
remove the first line of the text file
"""
lines = open('uncertainty_a.txt').readlines()
open('uncertainty_a_post.txt','w').writelines(lines[1:])
adata = np.loadtxt('uncertainty_a_post.txt')



"""
prepate arrays for the plots
"""
p = np.loadtxt('p.txt')
T = np.loadtxt('T.txt')
cm = plt.cm.get_cmap('plasma_r')
sc = plt.scatter(T, p,c=adata[:,2],vmin=adata[:,2].min(),\
                 vmax=adata[:,2].max(),s=30,cmap=cm, edgecolors='none')
plt.xlabel('Temperature [K]');plt.ylabel('Pressure [Pa]')
cbar = plt.colorbar(sc); cbar.set_label("Absolute Uncertainty")
plt.xlim([273.15,1073.15]);plt.ylim([0,100e6])
plt.grid(True)
plt.show()
