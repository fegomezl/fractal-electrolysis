import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec as gs
from matplotlib import cm

def load_data(config, ii):
    try:
        fields = np.loadtxt('results/data/data_'+str(ii)+'/fields_0.dat',delimiter=',')
        bit_map = np.loadtxt('results/data/data_'+str(ii)+'/bit_map_0.dat',delimiter=',')
        t = np.loadtxt('results/data/data_'+str(ii)+'/time.txt')
        for jj in range(1, config.nproc):
            fields = np.concatenate((fields, np.loadtxt('results/data/data_'+str(ii)+'/fields_'+str(jj)+'.dat',delimiter=',')))
            bit_map = np.concatenate((bit_map, np.loadtxt('results/data/data_'+str(ii)+'/bit_map_'+str(jj)+'.dat',delimiter=',')))
    except:
        try:
            ii = 0
            fields = np.loadtxt('results/data/data_'+str(ii)+'/fields_0.dat',delimiter=',')
            bit_map = np.loadtxt('results/data/data_'+str(ii)+'/bit_map_0.dat',delimiter=',')
            t = np.loadtxt('results/data/data_'+str(ii)+'/time.txt')
            for jj in range(1, nproc):
                fields = np.concatenate((fields, np.loadtxt('results/data/data_'+str(ii)+'/fields_'+str(jj)+'.dat',delimiter=',')))
                bit_map = np.concatenate((bit_map, np.loadtxt('results/data/data_'+str(ii)+'/fields_'+str(jj)+'.dat',delimiter=',')))
            print('Set default to initial data.')
        except:
            print('No data.')
            exit()

    return fields, bit_map, t

class Config(object):

    def __init__(self):
        #Parameters of the program
        parameters = open('settings/parameters.txt')

        for ii in range(0, 14):
            line = parameters.readline()

        line = parameters.readline()
        line = line.split("#")
        self.nproc = int(line[0])

        for ii in range(0, 1):
            line = parameters.readline()

        line = parameters.readline()
        line = line.split("#")
        self.iterations = int(line[0])

        line = parameters.readline()
        line = line.split("#")
        self.vis_iterations = int(line[0])

        self.frames = int(self.iterations/self.vis_iterations)

        line = parameters.readline()
        line = line.split("#")
        self.n = int(line[0])

        for ii in range(0, 8):
            line = parameters.readline()

        line = parameters.readline()
        line = line.split("#")
        self.dt = float(line[0])*self.vis_iterations

        line = parameters.readline()
        line = line.split("#")
        self.L = float(line[0])

        for ii in range(0, 18):
            line = parameters.readline()
        line = line.split("#")
        self.fps = int(line[0])

        parameters.close()
        
def plot_fileds_and_particles(config,fields,bit_map,X,Y,U,V,ii,d3,t):
    #Transform data
    domain = fields[:,0]
    phi = fields[:,1]
    electric_field_x = fields[:,2]
    electric_field_y = fields[:,3]
    density = fields[:,4]
    electric_field = np.hypot(electric_field_x, electric_field_y)

    Domain = np.reshape(domain, (-1, config.n))
    bit_map = np.reshape(bit_map, (-1, config.n))
    Phi = np.reshape(phi, (-1, config.n))
    Electric_field_x = np.reshape(electric_field_x, (-1, config.n))
    Electric_field_y = np.reshape(electric_field_y, (-1, config.n))
    density = np.reshape(density, (-1, config.n))
    Electric_field = np.reshape(electric_field, (-1, config.n))

    #Plot fields and particles
    fig = plt.figure()
    grid = gs.GridSpec(3, 2)

    fig.suptitle('t = '+str(t)+' s')

    #Particles and domain
    ax = plt.subplot(grid[:-1, :])
    ax.set(xlim=(-config.L/2, config.L/2), ylim=(-config.L/2, config.L/2))
    plt.pcolormesh(X, Y, density, cmap = cm.Blues)
    plt.colorbar()

    c_white = cm.colors.colorConverter.to_rgba('white',alpha = 0)
    c_black= cm.colors.colorConverter.to_rgba('black',alpha = 1)
    cmap_rb = cm.colors.LinearSegmentedColormap.from_list('rb_cmap',[c_black,c_white],255)

    plt.pcolormesh(X, Y, Domain, cmap=cmap_rb)

    #Electric potential
    ax = plt.subplot(grid[2, 0])
    ax.set(xlim=(-config.L/2, config.L/2), ylim=(-config.L/2, config.L/2))
    plt.pcolormesh(X, Y, Phi, cmap = cm.Blues)
    plt.colorbar()

    #Electric field
    ax = plt.subplot(grid[2, 1])
    ax.set(xlim=(-config.L/2, config.L/2), ylim=(-config.L/2, config.L/2))
    plt.streamplot(U, V, Electric_field_x, Electric_field_y, linewidth=0.5, arrowsize=0.5)
    plt.pcolormesh(X, Y, Electric_field, cmap = cm.Reds)
    plt.colorbar()

    plt.show()

    #Color map
    plt.pcolormesh(X, Y, bit_map, cmap = cm.binary_r)
    plt.show()

    if (d3 != '3d'):
        exit()

    #3D Electric potential
    fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    ax.set_xlim([-config.L/2, config.L/2])
    ax.set_ylim([-config.L/2, config.L/2])
    surf=ax.plot_surface(U, V, Phi, cmap=cm.Blues)
    cbar = fig.colorbar(surf)
    cbar.set_label('V')
    plt.show()
