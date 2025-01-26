
from calendar import c
from re import X
from struct import unpack
from this import d
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# u_R(r) 
rall = []
uall = []
rall, uall = np.loadtxt('ftrapezoidN4800.txt', usecols=(0, 1), unpack=True)
#plt.scatter(rall, uall, s=5, facecolors='none', edgecolors='red', label=r'$u_R(r)$')

# repulsive wca part
r, u0 = np.loadtxt('traprep.txt', usecols=(0, 1), unpack=True)
lam = 0.0575235
rep = u0 + lam 
#log = np.log(rep)
#print('For r=2.9 the pot value is:', log[201])
#plt.scatter(r, log, s=3, facecolors='none', edgecolors='blue', label=r'$\ln{(u_R^{repulsive}(r) / \epsilon_{GB})}$')
plt.scatter(r, rep, s=3, facecolors='none', edgecolors='blue', label=r'$u_R^{repulsive}(r) / \epsilon_{GB}$')
rem1 = np.linspace(2.951,3.7,100)
#plt.plot(rem1, 0*rem1, color='blue', lw=3)

# attractive wca part
rprime = []
uattr  = []
rprime, uattr = np.loadtxt('trapatt.txt', usecols=(0, 1), unpack=True)
#plt.scatter(rprime, uattr, s=3, facecolors='none', edgecolors='red', label=r'$u_R^{attractive}(r) / \epsilon_{GB}$')
rem2 = np.linspace(0,2.951,100)
#plt.plot(rem2, -lam + 0*rem2, color='red', lw=3)

# vertical line @r=2.95\sigma_{GB}
#plt.vlines(2.951, -14, 0, color='black', linestyle='--', lw=1, label=r'$r_{min}=2.951 \sigma_{GB}$')


# highlight -\lambda
#y_ticks = [-lam, -0.1, 0.6, 0, 0.1, 0.2, 0.3, 0.4, 0.5] 
#y_labels = [r'$-\lambda$', r'$-0.1$', 0.6, 0, 0.1, 0.2, 0.3, 0.4, 0.5]
#plt.yticks(y_ticks, y_labels)

# highlight r_{min}
#x_ticks = [2.951, 2.6, 2.7, 2.75, 2.85, 2.8, 2.9, 3, 3.2, 3.4, 3.7]
#x_labels = [r'$r_{min}$', 2.6, 2.7, 2.75, 2.85, 2.8, 2.9, 3, 3.2, 3.4, 3.7]
#plt.xticks(x_ticks, x_labels)


# Non-linear fit to the repulsive wall
rmin = 2.951
def fit(x, a, b):
   return a*np.exp(-(x/b)) - a*np.exp(-(rmin/b))
params, covs = curve_fit(fit, r, rep, maxfev=1000000)
print("params for original fit: ", params)
a, b = params[0], params[1]
yfit = a*np.exp(-(r/b)) - a*np.exp(-(rmin/b))
#lnfit = np.log(yfit)
#print('For r=2.9 the fit value is:', lnfit[201])
#plt.scatter(r, lnfit, s=3, facecolors='none', edgecolors='red', label=r'$\ln{(f(r)/\epsilon_{GB})}$')
plt.scatter(r, yfit, s=3, facecolors='none', edgecolors='red', label=r'$f(r)/\epsilon_{GB}$')

# Modified Non-linear fit to the repulsive wall
def modfit(x, a2, b2):
   return a2*np.exp(-((x-2.5)/b2)) - a2*np.exp(-((rmin-2.5)/b2))
params2, covs2 = curve_fit(modfit, r, rep, maxfev=1000000)
print("params for modified fit: ", params2)
a2, b2 = params2[0], params2[1]
yfit2 = a2*np.exp(-((r-2.5)/b2)) - a2*np.exp(-((rmin-2.5)/b2))
#lnfit = np.log(yfit)
#print('For r=2.9 the fit value is:', lnfit[201])
#plt.scatter(r, lnfit, s=3, facecolors='none', edgecolors='red', label=r'$\ln{(f(r)/\epsilon_{GB})}$')
plt.scatter(r, yfit2, s=3, facecolors='none', edgecolors='green', label=r'$f_2(r)/\epsilon_{GB}$')

# Plotting

plt.xlim(2.7,3)
plt.ylim(-0.1,1)

#plt.title(r'WCA decomposition of $u_R(r)$')
plt.title('Fit to the repulsive wall')
#plt.title('The corresponding logarithms')
xaxis = np.linspace(0,6,1000)
plt.plot(xaxis, 0*xaxis, color='black', lw=1)
plt.xlabel(r'$ r/\sigma_{GB} $', fontsize='20')
plt.tick_params(labelsize=15)
#plt.legend(loc='center left', fancybox=True, shadow=True)
plt.legend(loc='upper right', fancybox=True, shadow=True)
plt.show()