import sys
from tools import *

def main():
    #Parameters of the program
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
    ii = int(sys.argv[1])
    fields, particles = load_data(config,ii)
<<<<<<< HEAD

    #Plot data
    try:
        D3 = sys.argv[2]
    except:
        D3 = None

    plot_fileds_and_particles(config,fields,particles,X,Y,U,V,ii,D3)
=======

    #Plot data
    plot_fileds_and_particles(config,fields,particles,X,Y,U,V,ii,sys.argv[2])
>>>>>>> 0e7fa727dc3b8788e38dff084257d01ea9b80888

if __name__ == '__main__':
    main()