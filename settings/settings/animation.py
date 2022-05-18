import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.animation as animation
import matplotlib as mpl
from tools import *

def init():
    quad1.set_array([])
    scatter1.set_offsets([])
    return quad1,scatter1

def animate(n,config,quad1,scatter1):
    #Load data
    fields, particles = load_data(config,n)

    #Transform data
    domain = fields[:,0]
    Domain = np.reshape(domain, (-1, config.n))
    quad1.set_array(Domain.ravel())
    scatter1.set_offsets(particles[:,:2])

    return quad1,scatter1

def main():
    config = Config()
    frn = 10
    fps = 30

    #Grid
    x = np.linspace(-config.L/2, config.L/2, config.n+1)
    y = np.linspace(-config.L/2, config.L/2, config.n+1)
    X, Y = np.meshgrid(x, y)

    #Nodes
    u = np.linspace(-(config.L-config.L/config.n)/2, (config.L-config.L/config.n)/2, config.n)
    v = np.linspace(-(config.L-config.L/config.n)/2, (config.L-config.L/config.n)/2, config.n)
    U, V = np.meshgrid(u, v)

    #Load data
    ii = 0
    fields, particles = load_data(config,ii)

    #Transform data
    domain = fields[:,0]
    phi = fields[:,1]
    electric_field_x = fields[:,2]
    electric_field_y = fields[:,3]
    electric_field = np.hypot(electric_field_x, electric_field_y)

    Domain = np.reshape(domain, (-1, config.n))
    Phi = np.reshape(phi, (-1, config.n))
    Electric_field_x = np.reshape(electric_field_x, (-1, config.n))
    Electric_field_y = np.reshape(electric_field_y, (-1, config.n))
    Electric_field = np.reshape(electric_field, (-1, config.n))

    pos_x = np.append(particles[:,0], [2*config.L])
    pos_y = np.append(particles[:,1], [2*config.L])

    #create figure 
    fig = plt.figure(figsize=(10,10),facecolor='white')
    fig.suptitle('t = '+str(ii*config.vis_iterations*config.dt)+' s')

    #Particles and domain
    ax1 = plt.subplot()
    ax1.set(xlim=(-config.L/2, config.L/2), ylim=(-config.L/2, config.L/2))
    quad1 = plt.pcolormesh(X, Y, Domain, cmap = cm.binary_r) 
    scatter1 = plt.scatter(pos_x, pos_y, marker='o', s=1)

    anim = animation.FuncAnimation(fig,animate,fargs=(config,quad1,scatter1),frames=frn,interval=1,blit=False,repeat=False)
    anim.save('animation.gif', writer='imagemagick', fps=fps)
    #plt.show()

if __name__ == '__main__':
    main()
