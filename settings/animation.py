import matplotlib.pyplot as plt
import matplotlib.animation as animation

from tools import *

def init():
    quad1.set_array([])
    quad2.set_array([])
    return quad1,quad2

def animate(n,config,quad1,quad2,ax):
    #Load data
    fields, bit_map, t = load_data(config,n)
    ax.set_title('t = '+str(t)+' s')

    #Transform data
    Domain = fields[:,0]
    Domain = np.reshape(Domain, (-1, config.n))

    density = fields[:,4]
    density = np.reshape(density, (-1, config.n))

    #Update plot
    quad1.set_array(Domain.ravel())
    quad2.set_array(density.ravel())
    
    return quad1

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
    fields, bit_map, t = load_data(config,ii)

    #Transform data
    Domain = fields[:,0]
    Domain = np.reshape(Domain, (-1, config.n))

    density = fields[:,4]
    density = np.reshape(density, (-1, config.n))

    #Create figure 
    fig = plt.figure(facecolor='white')

    #Create transparency range 
    c_white = cm.colors.colorConverter.to_rgba('white',alpha = 0)
    c_black= cm.colors.colorConverter.to_rgba('black',alpha = 1)
    cmap_rb = cm.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_black,c_white],255)

    #Particles and domain
    ax1 = plt.subplot()
    ax1.set_title('t = '+str(t)+' s')
    ax1.set(xlim=(-config.L/2, config.L/2), ylim=(-config.L/2, config.L/2))
    quad2 = plt.pcolormesh(X, Y, density, cmap = cm.Blues) 
    clb = plt.colorbar()
    clb.set_label('Particle Density')
    quad1 = plt.pcolormesh(X, Y, Domain, cmap=cmap_rb)

    anim = animation.FuncAnimation(fig,animate,fargs=(config,quad1,quad2,ax1,),frames=26,interval=1,blit=False,repeat=False)
    anim.save('animation.gif', writer='imagemagick', fps=config.fps)


if __name__ == '__main__':
    main()
