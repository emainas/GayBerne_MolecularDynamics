import numpy as np
import matplotlib.pyplot as plt

r = []
c22 = []
c44 = []
c66 = []
c88 = []
c11 = []

rnew = []
trap = []

rnew, trap = np.loadtxt('ftrapezoidN4800.txt', usecols=(0, 1), unpack=True)

r, c22, c44, c66, c88, c11 = np.loadtxt('fcnn.txt', usecols=(0, 2, 3, 4, 5, 6), unpack=True)

def f(x):
    return x*0 + 0

z = np.linspace(0,6,100)

plt.plot(z, f(z), color='black', lw=1)

#plt.plot(rnew, trap, label=r'$u_R(r)/\epsilon_{GB}$', color='red', lw=2)

plt.plot(r, 4*c22, label=r'$4 c_{22}(r)/ \epsilon_{GB}$', color='red', lw=2)
plt.plot(r, 16*c44, label=r'$16 c_{44}(r)/ \epsilon_{GB}$', color='blue', lw=2)
plt.plot(r, 36*c66, label=r'$36 c_{66}(r)/ \epsilon_{GB}$', color='black', lw=2)
plt.plot(r, 64*c88, label=r'$64 c_{88}(r)/ \epsilon_{GB}$', color='cyan', lw=2)
plt.plot(r, 100*c11, label=r'$100 c_{1010}(r)/ \epsilon_{GB}$', color='green', lw=2)


plt.xlim(2.5,3.75)
plt.ylim(-0.1,3)

plt.xlabel(r'$ r/\sigma_{GB} $', fontsize='20')
plt.tick_params(labelsize=15)
plt.legend(loc='upper center')
plt.show()
