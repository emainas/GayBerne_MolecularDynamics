

import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('agg')

def side_to_side(x):
    
    epsilon_ss = 5/3
    sigma_ss = 0

    return 4 * epsilon_ss * ( ( 1 / (x+sigma_ss) )**12 - ( 1 / (x+sigma_ss))**6)

def side_to_end(x):
    
    chi = 0.8
    chi_prime = 0.3819660112501052084965635913249570876359939575195312

    epsilon_se = (1 - chi_prime)**2
    sigma_se = 1 - (1 - chi)**(-1/2)

    return 4 * epsilon_se * ( ( 1 / (x+sigma_se) )**12 - ( 1 / (x+sigma_se))**6)

def end_to_end(x):
    
    chi = 0.8
    chi_prime = 0.3819660112501052084965635913249570876359939575195312
    
    epsilon_ee = (1/np.sqrt(1 - chi**2)) * (1 - (chi_prime/2)*(4/(1 + chi_prime)) )**2
    sigma_ee = 1 - (1 - (chi/2)*(4/(1 + chi)))**(-1/2)

    return 4 * epsilon_ee * ((1/(x+sigma_ee))**12 - (1/(x+sigma_ee))**6)

x_se = np.linspace(1, 5, 10000)
y_se = side_to_end(x_se)

x_ss = np.linspace(0.1, 5, 10000)
y_ss = side_to_side(x_ss)

x_ee = np.linspace(2.8, 5, 10000)
y_ee = end_to_end(x_ee)

plt.plot(x_ss, y_ss, color='blue', label='Side to side', lw=4)
#plt.plot(x_se, y_se, color='green', label='Side to end', lw=4)
plt.plot(x_ee, y_ee, color='red', label='End to end', lw=4)

min_ss_index = np.argmin(y_ss)
min_ss_x = x_ss[min_ss_index]
min_ss_y = y_ss[min_ss_index]

min_ee_index = np.argmin(y_ee)
min_ee_x = x_ee[min_ee_index]
min_ee_y = y_ee[min_ee_index]

plt.hlines(y=min_ss_y, xmin=0, xmax=min_ss_x, linestyle='dashed', color='blue')
plt.hlines(y=min_ee_y, xmin=0, xmax=min_ee_x, linestyle='dashed', color='red')

plt.vlines(x=1, ymin=0, ymax=-2, linestyle='dashed', color='blue')
plt.vlines(x=3, ymin=0, ymax=-2, linestyle='dashed', color='red')

y_ticks = [-5*min_ee_y, -3*min_ee_y, -min_ee_y, 0, min_ee_y, 3*min_ee_y, 5*min_ee_y]
y_labels = [r'$5/3$', r'$1$', r'$1/3$', r'$0$', r'$-1/3$', r'$-1$', r'$-5/3$']
plt.yticks(y_ticks, y_labels)

ax = plt.gca()  # Get the current axes
line_width = 3
ax.spines['top'].set_linewidth(line_width)
ax.spines['bottom'].set_linewidth(line_width)
ax.spines['left'].set_linewidth(line_width)
ax.spines['right'].set_linewidth(line_width)

plt.axhline(y=0, xmin=0, xmax=5, lw=2, color='black')
plt.xlim(0, 5)
plt.ylim(-2, 2)
plt.xlabel(r'$r / \sigma_0$', fontsize='30')
plt.ylabel(r'$U / \epsilon_0$', fontsize='25') 
plt.tick_params(labelsize=20)
#plt.legend(loc='upper right', prop={'size': 18}, frameon=False) 
plt.savefig(fr'gay_berne.jpeg', dpi=600, bbox_inches='tight')