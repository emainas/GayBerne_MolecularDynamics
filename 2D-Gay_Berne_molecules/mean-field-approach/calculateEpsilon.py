
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

r = []
c22 = []
c44 = []
c66 = []
c88 = []
c1010 = []
g = []

r, c22, c44, c66, c88, c1010 = np.loadtxt('fcnn.txt', usecols=(0, 2, 3, 4, 5, 6), unpack=True)
g = np.loadtxt('g.txt', usecols=(1), unpack=True)

prod22 = 4*c22*g*r
prod44 = 16*c44*g*r
prod66 = 36*c66*g*r
prod88 = 64*c88*g*r
prod1010 = 100*c1010*g*r

plt.plot(r, prod22, color='red', label=r'$r 2^2 g(r) (c_{22}(r)/ \epsilon_{GB})$', lw=2)
plt.plot(r, prod44, color='blue', label=r'$r 4^2 g(r) (c_{44}(r)/ \epsilon_{GB})$', lw=2)
plt.plot(r, prod66, color='black', label=r'$r 6^2 g(r) (c_{66}(r)/ \epsilon_{GB})$', lw=2)
#plt.plot(r, prod88, color='green', label=r'$8^2 g(r) c_{88}(r)$')
#plt.plot(r, prod1010, color='cyan', label=r'$10^2 g(r) c_{1010}(r)$')

#Integrating the products

int22 = 0
int44 = 0
int66 = 0
int88 = 0
int1010 = 0
h = 0.001
rhoR = 0.035648
for i in range(3999):
    int22 += (1/2)*h*(prod22[i] + prod22[i+1])
    int44 += (1/2)*h*(prod44[i] + prod44[i+1])
    int66 += (1/2)*h*(prod66[i] + prod66[i+1])
    int88 += (1/2)*h*(prod88[i] + prod88[i+1])
    int1010 += (1/2)*h*(prod1010[i] + prod1010[i+1])

epsilon = (-1/8)*rhoR*2*np.pi*(int22 + int44 + int66)
print('Trapezoidal rule epsilon is:', epsilon)

#Simpson rule
temp = integrate.simps(prod22 + prod44 + prod66, r)
epsilonsimp = (-1/8)*2*np.pi*rhoR*temp
print('Simpson rule epsilon is', epsilonsimp)

xaxis = np.linspace(2,4,1000)
plt.plot(xaxis, 0*xaxis, lw=1, color='black')
plt.xlim(2,4)
plt.ylim(-0.5,2.5)
plt.xlabel(r'$r/ \sigma_{GB}$', fontsize='20')
plt.legend(loc='upper right')
plt.tick_params(labelsize=15)
plt.show()