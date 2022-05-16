import sys
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec as gs
from matplotlib import cm

#Parameters of the program
parameters = open('settings/parameters.txt')

for ii in range(0, 14):
    line = parameters.readline()

line = parameters.readline()
line = line.split("#")
nproc = int(line[0])

for ii in range(0, 2):
    line = parameters.readline()

line = parameters.readline()
line = line.split("#")
vis_iterations = int(line[0])

line = parameters.readline()
line = line.split("#")
n = int(line[0])

for ii in range(0, 8):
    line = parameters.readline()

line = parameters.readline()
line = line.split("#")
dt = float(line[0])*vis_iterations

line = parameters.readline()
line = line.split("#")
L = float(line[0])

parameters.close()

#Grid
x = np.linspace(-L/2, L/2, n+1)
y = np.linspace(-L/2, L/2, n+1)
X, Y = np.meshgrid(x, y)

#Nodes
u = np.linspace(-(L-L/n)/2, (L-L/n)/2, n)
v = np.linspace(-(L-L/n)/2, (L-L/n)/2, n)
U, V = np.meshgrid(u, v)

#Load data
try:
    ii = int(sys.argv[1])
    fields = np.loadtxt('results/data/data_'+str(ii)+'/fields_0.dat')
    particles = np.loadtxt('results/data/data_'+str(ii)+'/particles_0.dat')
    for jj in range(1, nproc):
        fields = np.concatenate((fields, np.loadtxt('results/data/data_'+str(ii)+'/fields_'+str(jj)+'.dat')))
        particles = np.concatenate((particles, np.loadtxt('results/data/data_'+str(ii)+'/particles_'+str(jj)+'.dat')))
except:
    try:
        ii = 0
        fields = np.loadtxt('results/data/data_'+str(ii)+'/fields_0.dat')
        particles = np.loadtxt('results/data/data_'+str(ii)+'/particles_0.dat')
        for jj in range(1, nproc):
            fields = np.concatenate((fields, np.loadtxt('results/data/data_'+str(ii)+'/fields_'+str(jj)+'.dat')))
            particles = np.concatenate((particles, np.loadtxt('results/data/data_'+str(ii)+'/particles_'+str(jj)+'.dat')))
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

Domain = np.reshape(domain, (-1, n))
Phi = np.reshape(phi, (-1, n))
Electric_field_x = np.reshape(electric_field_x, (-1, n))
Electric_field_y = np.reshape(electric_field_y, (-1, n))
Electric_field = np.reshape(electric_field, (-1, n))

try:
    pos_x = np.append(particles[:,0], [2*L])
    pos_y = np.append(particles[:,1], [2*L])
except:
    pos_x = [2*L]
    pos_y = [2*L]

#Plot fields and particles
fig = plt.figure(figsize=(10,10))
grid = gs.GridSpec(3, 2)

fig.suptitle('t = '+str(ii*dt)+' s')

#Particles and domain
ax = plt.subplot(grid[:-1, :])
ax.set(xlim=(-L/2, L/2), ylim=(-L/2, L/2))
plt.pcolormesh(X, Y, Domain, cmap = cm.binary_r)
plt.scatter(pos_x, pos_y, marker='o', s=1)

#Electric potential
ax = plt.subplot(grid[2, 0])
ax.set(xlim=(-L/2, L/2), ylim=(-L/2, L/2))
plt.pcolormesh(X, Y, Phi, cmap = cm.Blues)
plt.colorbar()

#Electric field
ax = plt.subplot(grid[2, 1])
ax.set(xlim=(-L/2, L/2), ylim=(-L/2, L/2))
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
ax.set_xlim([-L/2, L/2])
ax.set_ylim([-L/2, L/2])
surf=ax.plot_surface(U, V, Phi, cmap=cm.Blues)
cbar = fig.colorbar(surf)
cbar.set_label('V')
plt.show()
