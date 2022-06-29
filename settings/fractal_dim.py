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

    fields, particles, bit_map, t = load_data(config,n)
    bit_map = 1 - np.reshape(bit_map, (-1, config.n))
    image = np.array(bit_map * 255, dtype = np.uint8)

    #plt.imshow(image)
    #plt.show()

    data = ps.metrics.boxcount(image,bins=np.array([6,8,12,16,20,25,32,44,50]))

    L = np.log10(data.size)
    N = np.log10(data.count)

    result = stats.linregress(L,N)
    print('Fractial dimention = {:.4f} +- {:.4f}'.format(-result.slope, result.stderr))
    
    '''
    label = r'$log(N)= D log(L)+b$''\n'r'$D$''={:.3f}±{:.3f} \n b={:.3f}±{:.3f} \n r={:.5f}'.format(result.slope, result.stderr, result.intercept, result.intercept_stderr, result.rvalue)
    fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))
    ax1.set_xlabel('log(box size)')
    ax1.set_ylabel('log(count)')
    ax1.plot(L, N,'-o')
    ax1.plot(L,result.slope*L+result.intercept,color='red',label=label)
    ax1.legend()

    plt.show()
    '''

if __name__ == '__main__':
    main()
