import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm

#Parameters of the program
parameters = open('settings/parameters.txt')

for ii in range(0, 16):
    line = parameters.readline()

line = parameters.readline()
line = line.split("#")
vis_iterations = int(line[0])

line = parameters.readline()
line = line.split("#")
nx = ny = int(line[0])

for ii in range(0, 8):
    line = parameters.readline()

line = parameters.readline()
line = line.split("#")
dt = float(line[0])*vis_iterations

line = parameters.readline()
line = line.split("#")
Lx = Ly = float(line[0])

parameters.close()

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

try:
    pos_x = np.append(particles[:,0], [2*Lx])
    pos_y = np.append(particles[:,1], [2*Ly])
except:
    pos_x = [2*Lx]
    pos_y = [2*Ly]

#Plot fields and particles
fig = plt.figure(figsize=(10,10))
gs = gridspec.GridSpec(3, 2)

fig.suptitle('t = '+str(n*dt)+' s')

#Particles and domain
ax = plt.subplot(gs[:-1, :])
ax.set(xlim=(-Lx/2, Lx/2), ylim=(-Lx/2, Lx/2))
plt.pcolormesh(X, Y, Domain, cmap = cm.binary_r)
plt.scatter(pos_x, pos_y, marker='o', s=1)

#Electric potential
ax = plt.subplot(gs[2, 0])
ax.set(xlim=(-Lx/2, Lx/2), ylim=(-Lx/2, Lx/2))
plt.pcolormesh(X, Y, Phi, cmap = cm.Blues)
plt.colorbar()

#Electric field
ax = plt.subplot(gs[2, 1])
ax.set(xlim=(-Lx/2, Lx/2), ylim=(-Lx/2, Lx/2))
plt.streamplot(U, V, Electric_field_x, Electric_field_y, linewidth=0.5, arrowsize=0.5)
plt.pcolormesh(X, Y, Electric_field, cmap = cm.Reds)
plt.colorbar()

plt.show()

try:
    D3 = sys.argv[2]
    if (D3 != '3d'):
        exit()
except:
    exit()

#3D Electric potential
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
ax.set_xlim([-Lx/2, Lx/2])
ax.set_ylim([-Ly/2, Ly/2])
surf=ax.plot_surface(U, V, Phi, cmap=cm.Blues)
cbar = fig.colorbar(surf)
cbar.set_label('V')
plt.show()
