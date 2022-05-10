import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#Parameters from the domain
Lx = 16.
Ly = 16.
nx = 729
ny = 729

#Grid
x = np.linspace(-Lx/2, Lx/2, nx+1) 
y = np.linspace(-Ly/2, Ly/2, ny+1) 
X, Y = np.meshgrid(x, y)

#Nodes
u = np.linspace(-(Lx-Lx/nx)/2, (Lx-Lx/nx)/2, nx) 
v = np.linspace(-(Ly-Ly/ny)/2, (Ly-Ly/ny)/2, ny) 
U, V = np.meshgrid(u, v)

#Particles and boundary
z = np.loadtxt('results/boundary.dat')
Z = np.reshape(z, (-1, nx))

z = np.loadtxt('results/particles.dat')
zx = z[:,0]
zy = z[:,1]

plt.title('Particles and Electrodes')
plt.pcolormesh(X, Y, Z, cmap = cm.binary_r) 
plt.scatter(zx, zy, marker='o', s=1)
plt.show()

#Electric potential
z = np.loadtxt('results/phi.dat')
Z = np.reshape(z, (-1, nx))

#2D
plt.title('Electric potential')
plt.pcolormesh(X, Y, Z, cmap = cm.Blues) 
cbar = plt.colorbar()
cbar.set_label('V'); 
plt.show()

#3D
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set_xlim([-Lx/2, Lx/2])
ax.set_ylim([-Ly/2, Ly/2])
surf=ax.plot_surface(U, V, Z, cmap=cm.Blues)
cbar = fig.colorbar(surf)
cbar.set_label('V'); 
plt.show()

#Electric field
zx = np.loadtxt('results/gradientx.dat')
zy = np.loadtxt('results/gradienty.dat')
Zx = np.reshape(zx, (-1, nx))
Zy = np.reshape(zy, (-1, nx))

z = np.hypot(zx, zy)
Z = np.reshape(z, (-1, nx))

plt.title('Electric field')
plt.streamplot(U, V, Zx, Zy)
plt.pcolormesh(X, Y, Z, cmap = cm.Reds)
cbar = plt.colorbar()
cbar.set_label('V/cm'); 
plt.show()
