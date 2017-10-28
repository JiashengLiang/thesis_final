import numpy as np
import matplotlib.pyplot as plt

lines = open('uncertainty_p.txt').readlines()
open('uncertainty_p_post.txt','w').writelines(lines[1:])
adata = np.loadtxt('uncertainty_p_post.txt')
p = np.loadtxt('p.txt')
T = np.loadtxt('T.txt')
cm = plt.cm.get_cmap('plasma_r')
sc = plt.scatter(T, p,c=adata[:,3],vmin=adata[:,3].min(),\
                 vmax=adata[:,3].max(),s=30,cmap=cm, edgecolors='none')
plt.xlabel('Temperature [K]');plt.ylabel('Pressure [Pa]')
plt.title('Pressure Calculated from updateThermoFromRHOU()',y=1.05)
cbar = plt.colorbar(sc); cbar.set_label("Percentage Uncertainty [%]")
plt.xlim([273.15,1073.15]);plt.ylim([0,100e6])
plt.grid(True)
plt.show()
