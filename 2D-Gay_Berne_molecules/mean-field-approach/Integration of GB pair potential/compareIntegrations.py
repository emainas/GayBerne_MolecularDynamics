
from re import X
from struct import unpack
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

r = []
trap = []
compsimp = []
boole = []

r2 = []
trap2 = []
compsimp2 = []
boole2 = []

r, trap = np.loadtxt('ftrapezoidN20.txt', usecols=(0, 1), unpack=True)
compsimp = np.loadtxt('fcompsimpN20.txt', usecols=(1), unpack=True)
boole = np.loadtxt('fbooleN20.txt', usecols=(1), unpack=True)

r2, trap2 = np.loadtxt('ftrapezoidN4800.txt', usecols=(0, 1), unpack=True)
compsimp2 = np.loadtxt('fcompsimpN4800.txt', usecols=(1), unpack=True)
boole2 = np.loadtxt('fbooleN4800.txt', usecols=(1), unpack=True)

z = np.linspace(0,6,1000)

plt.plot(z, 0*z, color='black', lw=1)

#plt.scatter(r, trap, s=20, marker='x', color='red', label='Trapezoid rule, n=20')
#plt.scatter(r, compsimp, s=10, marker='x', color='blue', label='Composite Simpson rule, n=20')
#plt.scatter(r, boole, s=1, marker='x', color='black', label='Boole rule, n=20')

plt.scatter(r2, trap2, s=3, marker='o', color ='red', label=r'Trapezoid rule, $200 \times 200$')
#plt.scatter(r2, compsimp2, s=10, marker='x', color ='blue', label='Simpson rule, n=4800')
#plt.scatter(r2, boole2, s=1, marker='x', color ='black', label='Boole rule, n=4800')

# Customize y-axis and draw horizontal line
lam = 0.0575235
y_ticks = [-lam, -0.1, 0.6, 0, 0.1, 0.2, 0.3, 0.4, 0.5] 
y_labels = [r'$-\lambda$', r'$-0.1$', 0.6, 0, 0.1, 0.2, 0.3, 0.4, 0.5]
plt.yticks(y_ticks, y_labels)
temp = np.linspace(2.6, 2.951, 100)
plt.plot(temp, temp*0 - lam, color='black', linestyle='dotted', lw=2, label=r'$\lambda = 0.0575235 \epsilon_{GB}$')

x_ticks = [2.951, 2.8275, 2.6, 3, 3.2, 3.4, 3.6, 3.8]
x_labels = [r'$r_{min}$', r'$r_{cr}$', 2.6, 3, 3.2, 3.4, 3.6, 3.8]
plt.xticks(x_ticks, x_labels)

# vertical line @r=2.95\sigma_{GB}
plt.vlines(2.8275, -0.1, 0, color='black', linestyle='--', lw=1, label=r'$r_{cr}=2.8275 \sigma_{GB}$')
plt.vlines(2.951, -0.1, 0, color='black', linestyle='--', lw=1, label=r'$r_{min}=2.951 \sigma_{GB}$')



plt.xlim(2.6,4)
plt.ylim(-0.1,0.6)

plt.xlabel(r'$ r/\sigma_{GB} $', fontsize='20')
plt.ylabel(r'$  u_R(r)/\epsilon_{GB} $', fontsize='20') 
plt.tick_params(labelsize=15)
plt.legend(loc='upper right', fancybox=True, shadow=True)
plt.show()