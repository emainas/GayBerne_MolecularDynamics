
import matplotlib.pyplot as plt 
import numpy as np

x = [] 
y = [] 
s = []  
c = [] 



x, y = np.loadtxt('GB_pos.txt', usecols=(0, 1), unpack=True)

plt.scatter(x,y)

plt.show()


