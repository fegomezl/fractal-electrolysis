from tools import *
import os
import matplotlib.pyplot as plt
import porespy as ps
from scipy import stats
import numpy as np

def main():
    #Parameters of the program
    config = Config()

    #Load last data point
    path = 'results/data'
    n=len(next(os.walk(path))[1])-1

    fields, particles, bit_map = load_data(config,n)
    bit_map = 1 - np.reshape(bit_map, (-1, config.n))

    data = ps.metrics.boxcount(bit_map)

    result = stats.linregress(np.log10(data.size),np.log10(data.count))
    print('Fractial dimention = {:.4f} +- {:.4f}'.format(result.slope, result.stderr))
    label = r'$N=10^b L^d$''\n'r'$d$''={:.3f}±{:.3f} \n b={:.3f}±{:.3f} \n r={:.5f}'.format(result.slope, result.stderr, result.intercept, result.intercept_stderr, result.rvalue)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 4))
    ax1.set_yscale('log')
    ax1.set_xscale('log')
    ax1.set_xlabel('box edge length')
    ax1.set_ylabel('number of boxes spanning phases')
    ax2.set_xlabel('box edge length')
    ax2.set_ylabel('slope')
    ax2.set_xscale('log')
    ax1.plot(data.size, data.count,'-o')
    ax1.plot(data.size,np.power(10,result.intercept)*np.power(data.size,result.slope),color='red',label=label)
    ax2.plot(data.size, data.slope,'-o')
    ax1.legend()
    plt.show()

if __name__ == '__main__':
    main()