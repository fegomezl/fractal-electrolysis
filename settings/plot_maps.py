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
ax.plot_surface(x_index.reshape(ny,nx),y_index.reshape(ny,nx),z.reshape(ny,nx))
plt.title('Potential')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.imshow(z.reshape(ny,nx))
plt.title('Potential')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.imshow(zz.reshape(ny,nx))
plt.title('Boundary')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

plt.imshow(zzz.reshape(ny,nx))
plt.title('Dissociation')
plt.xlabel('x')
plt.ylabel('y')
plt.show()

'''
f = np.gradient(z.reshape(ny,nx))

plt.quiver(x_index.reshape(ny,nx),y_index.reshape(ny,nx),f[1],f[0])
plt.title('Numpy Gradient')
plt.xlabel('x')
plt.ylabel('y')
plt.show()
'''

plt.quiver(x_index.reshape(ny,nx),y_index.reshape(ny,nx),gradx.reshape(ny,nx),grady.reshape(ny,nx))
plt.title('Gradient')
plt.xlabel('x')
plt.ylabel('y')
plt.show()