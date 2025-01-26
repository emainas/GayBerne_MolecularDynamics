
from calendar import c
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

r = []
u0 = []

r, u0 = np.loadtxt('alldata.txt', usecols=(0, 1), unpack=True)

def rep(x, a, b, c):
   return a*np.exp(-b*x**3) - c*(1/x**6) 

params, covs = curve_fit(rep, r, u0, maxfev=1000000)
print("params for rep: ", params)
print("covariance for rep: ", covs)

a, b, c = params[0], params[1], params[2]
z = np.linspace(0,2.955,1000)
zprime = np.linspace(0,6,1000)

# a = 21806432.9 | b = 0.865031401 | c = 34.2271197

yrep = a*np.exp(-b*(r**3))  - c*(1/r**6) 

plt.plot(r, yrep, label=r'$fit$',color='blue', lw=1)
plt.plot(zprime, 0*z, color='black', lw=1)
plt.scatter(r, u0, s=1, facecolors='none', edgecolors='red', label='Reference potential energy data')

plt.xlim(2.6,4.2)
plt.ylim(-0.2,3)

plt.xlabel(r'$ r/\sigma_0 $', fontsize='20')
plt.ylabel(r'$  u_0(r)/\epsilon_0 $', fontsize='20') 
plt.tick_params(labelsize=15)
plt.legend(loc='upper right', fancybox=True, shadow=True)
plt.show()