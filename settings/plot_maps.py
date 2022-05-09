import numpy as np
import matplotlib.pyplot as plt

data=np.loadtxt('results/phi.txt')
x_index = data[:,0]
y_index = data[:,1]
z = data[:,2]

data=np.loadtxt('results/boundary.txt')
zz = data[:,2]


data=np.loadtxt('results/dissociation.txt')
zzz = data[:,2]

gradx=np.loadtxt('results/gradientx.txt')[:,2]
grady=np.loadtxt('results/gradienty.txt')[:,2]

nx = np.shape(np.unique(x_index))[0]
ny = np.shape(np.unique(y_index))[0]


ax = plt.axes(projection ='3d')
ax.plot_surface(x_index.reshape(nx,ny),y_index.reshape(nx,ny),z.reshape(nx,ny))
plt.title('Potential')
plt.show()

plt.imshow(z.reshape(nx,ny))
plt.title('Potential')
plt.show()

plt.imshow(zz.reshape(nx,ny))
plt.title('Boundary')
plt.show()

plt.imshow(zzz.reshape(nx,ny))
plt.title('Dissociation')
plt.show()

plt.quiver(x_index,y_index,gradx,grady)
plt.title('Gradient')
plt.show()