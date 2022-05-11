import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#Parameters from the domain
Lx = 10.
Ly = 10.
nx = 729
ny = 729
dt = 1.

#Grid
x = np.linspace(-Lx/2, Lx/2, nx+1) 
y = np.linspace(-Ly/2, Ly/2, ny+1) 
X, Y = np.meshgrid(x, y)

#Nodes
u = np.linspace(-(Lx-Lx/nx)/2, (Lx-Lx/nx)/2, nx) 
v = np.linspace(-(Ly-Ly/ny)/2, (Ly-Ly/ny)/2, ny) 
U, V = np.meshgrid(u, v)

#Load data
try:
    n = int(sys.argv[1])
    fields = np.loadtxt('results/data/fields_'+str(n)+'.dat')
    particles = np.loadtxt('results/data/particles_'+str(n)+'.dat')
except:
    try:
        n = 0
        fields = np.loadtxt('results/data/fields_'+str(n)+'.dat')
        particles = np.loadtxt('results/data/particles_'+str(n)+'.dat')
        print('Set default to initial data.')
    except:
        print('No data.')
        exit()

#Transform data
domain = fields[:,0]
phi = fields[:,1]
electric_field_x = fields[:,2]
electric_field_y = fields[:,3]
electric_field = np.hypot(electric_field_x, electric_field_y)

Domain = np.reshape(domain, (-1, nx))
Phi = np.reshape(phi, (-1, nx))
Electric_field_x = np.reshape(electric_field_x, (-1, nx))
Electric_field_y = np.reshape(electric_field_y, (-1, nx))
Electric_field = np.reshape(electric_field, (-1, nx))

pos_x = particles[:,0]
pos_y = particles[:,1]

#Particles and boundary
plt.title('Particles and Electrodes (t = '+str(n*dt)+' s)')
plt.pcolormesh(X, Y, Domain, cmap = cm.binary_r) 
plt.scatter(pos_x, pos_y, marker='o', s=1)
plt.show()

#Electric potential
#2D
plt.title('Electric potential (t = '+str(n*dt)+' s)')
plt.pcolormesh(X, Y, Phi, cmap = cm.Blues) 
cbar = plt.colorbar()
cbar.set_label('V'); 
plt.show()
#3D
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set_xlim([-Lx/2, Lx/2])
ax.set_ylim([-Ly/2, Ly/2])
surf=ax.plot_surface(U, V, Phi, cmap=cm.Blues)
cbar = fig.colorbar(surf)
cbar.set_label('V'); 
plt.show()

#Electric field
plt.title('Electric field (t = '+str(n*dt)+' s)')
plt.streamplot(U, V, Electric_field_x, Electric_field_y)
plt.pcolormesh(X, Y, Electric_field, cmap = cm.Reds)
cbar = plt.colorbar()
cbar.set_label('V/cm'); 
plt.show()
