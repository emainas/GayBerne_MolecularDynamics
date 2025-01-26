
import matplotlib.pyplot as plt 
import numpy as np

x = []
y = []
z = []

x, y, w = np.loadtxt('energy1.txt', usecols=(0, 1, 2), unpack=True)

q = x/100

plt.plot(q, y/128, label="Potential", color='red')
plt.plot(q, w, label="Kinetic", color='blue')

#plt.xlim(0,350)
#plt.ylim(0,50)

plt.xlabel(r'$t / \tau$', fontsize='20')
plt.ylabel(r'$E/N \epsilon_0$', fontsize='20')
plt.tick_params(labelsize=15)
plt.legend(loc='upper right')
plt.show()



