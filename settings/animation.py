import matplotlib.pyplot as plt
import matplotlib.animation as animation

from tools import *

def init():
    quad1.set_array([])
    scatter1.set_offsets([])
    return quad1,scatter1

def animate(n,config,quad1,scatter1,ax):
    #Load data
    ax.set_title('t = '+str(n*config.dt)+' s')
    fields, particles, bit_map = load_data(config,n)

    #Transform data
    domain = fields[:,0]
    Domain = np.reshape(domain, (-1, config.n))

    #Update plot
    quad1.set_array(Domain.ravel())
    scatter1.set_offsets(particles[:,:2])

    return quad1,scatter1

def main():
    #Parameters from the domain
    config = Config()
    
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
    Domain = np.reshape(domain, (-1, config.n))

    try:
        pos_x = np.append(particles[:,0], [2*config.L])
        pos_y = np.append(particles[:,1], [2*config.L])
    except:
        pos_x = [2*config.L]
        pos_y = [2*config.L]

    #create figure 
    fig = plt.figure(figsize=(10,10),facecolor='white')

    #Particles and domain
    ax1 = plt.subplot()
    ax1.set_title('t = '+str(ii*config.dt)+' s')
    ax1.set(xlim=(-config.L/2, config.L/2), ylim=(-config.L/2, config.L/2))
    quad1 = plt.pcolormesh(X, Y, Domain, cmap = cm.binary_r) 
    scatter1 = plt.scatter(pos_x, pos_y, marker='o', s=1)

    anim = animation.FuncAnimation(fig,animate,fargs=(config,quad1,scatter1,ax1,),frames=config.frames,interval=1,blit=False,repeat=False)
    anim.save('animation.gif', writer='imagemagick', fps=config.fps)


if __name__ == '__main__':
    main()
