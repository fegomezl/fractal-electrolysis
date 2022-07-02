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
    fields, bit_map, t = load_data(config,ii)

    #Plot data
    try:
        d3=sys.argv[2]
    except:
        d3 = None

    plot_fileds_and_particles(config,fields,bit_map,X,Y,U,V,ii,d3,t)

if __name__ == '__main__':
    main()
