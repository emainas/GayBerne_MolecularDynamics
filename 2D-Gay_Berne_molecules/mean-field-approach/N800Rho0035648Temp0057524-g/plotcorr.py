
import matplotlib.pyplot as plt 
import numpy as np
from scipy import integrate

r_R = []
g_R = []
r_GB = []
g_GB = []
rho_R = 0.035648
rho_GB = 0.285

r_R, g_R = np.loadtxt('g.txt', usecols=(0, 1), unpack=True)
rscaled = r_R*(1/2.8275)
r_GB, g_GB = np.loadtxt('gGB.txt', usecols=(0, 1), unpack=True)

x_R, int_g_R = np.loadtxt('g_R_int.txt', usecols=(0, 1), unpack=True)
xscaled = x_R*(1/2.8275)
x_GB, int_g_GB = np.loadtxt('g_GB_int.txt', usecols=(0, 1), unpack=True)

int_R = integrate.simps(int_g_R, xscaled)
int_GB = integrate.simps(int_g_GB, x_GB)

print('The ref system coordination is', int_R)
print('The GB system coordination is', int_GB)


plt.plot(rscaled, g_R, label=r"$N=800, \rho_R^*=0.035648, T_R^*=0.057524$", color='red', lw=2)
plt.plot(r_GB, g_GB, label=r"$N=800, \rho_{GB}^*=0.285, T_{GB}^*=1$", color='blue', lw=2)

one = np.full_like(r_GB, 1)

plt.plot(r_GB, one, color='black', lw=1)
plt.xlim(0,5)
plt.ylim(0,1.5)
plt.xticks(np.arange(0,5,0.5))
plt.ylabel(r'$g(r)$', fontsize='20')
plt.xlabel(r'$r/ \sigma_{eff}$', fontsize='20')
plt.tick_params(labelsize=15)
plt.legend(loc='upper right')
plt.show()



