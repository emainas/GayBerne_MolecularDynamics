
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

N = 500 # number of particles
M = 1000 # number of frames
D = 3 # number of dimensions
Lx = Ly = Lz = 6.5

with open('trajort.bin', 'rb') as f:
    orientation_data = np.fromfile(f, dtype=np.float64) 

with open('trajpos.bin', 'rb') as f:
    position_data = np.fromfile(f, dtype=np.float64) 

def draw_vector_field(frame):
    
    ax.clear()
    ax.set_xlim(-Lx, Lx)
    ax.set_ylim(-Ly, Ly)
    ax.set_zlim(-Lz, Lz)

    x = np.zeros([N])
    y = np.zeros([N])
    z = np.zeros([N])
    u = np.zeros([N])
    v = np.zeros([N])
    w = np.zeros([N])

    for i in range(N):
        index = (i * D) + (frame * D * N) 
        print('index is:', index)
        print('i is:', i)
        x[i] = position_data[index]
        y[i] = position_data[index + 1]
        z[i] = position_data[index + 2]
        u[i] = orientation_data[index]
        v[i] = orientation_data[index + 1]
        w[i] = orientation_data[index + 2]

    ax.quiver(x, y, z, u, v, w, length=0.5, normalize=True, arrow_length_ratio=0, color='k')

    # plot black sticks over the quiver plot
    for i in range(N):
        ax.plot([x[i], x[i] + u[i]], [y[i], y[i] + v[i]], [z[i], z[i] + w[i]], color='black')


    ax.set_title('Frame {}'.format(frame))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ani = animation.FuncAnimation(fig, draw_vector_field, frames=range(M), interval=40)
ani.save('vector_field_animation.mp4', writer=animation.FFMpegWriter(fps=10))