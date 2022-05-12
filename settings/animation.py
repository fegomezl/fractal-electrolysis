import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
import matplotlib.animation as animation
import matplotlib as mpl

#Parameters from the domain
Lx = 10.
Ly = 10.
nx = 729
ny = 729
dt = 1.
vis_iterations = 60
dt = vis_iterations*dt
frn = 60
fps = 10
interval = 10


#Grid
x = np.linspace(-Lx/2, Lx/2, nx+1) 
y = np.linspace(-Ly/2, Ly/2, ny+1) 
X, Y = np.meshgrid(x, y)

#Nodes
u = np.linspace(-(Lx-Lx/nx)/2, (Lx-Lx/nx)/2, nx) 
v = np.linspace(-(Ly-Ly/ny)/2, (Ly-Ly/ny)/2, ny) 
U, V = np.meshgrid(u, v)

#Load data
n = 0
fields = np.loadtxt('results/data/fields_'+str(n)+'.dat')
particles = np.loadtxt('results/data/particles_'+str(n)+'.dat')

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

pos_x = np.append(particles[:,0], [2*Lx])
pos_y = np.append(particles[:,1], [2*Ly])

#create figure 
fig = plt.figure(figsize=(10,10),facecolor='white')
fig.suptitle('t = '+str(n*dt)+' s')
gs = gridspec.GridSpec(3, 2)

#Particles and domain
ax1 = plt.subplot(gs[:-1, :])
ax1.set(xlim=(-Lx/2, Lx/2), ylim=(-Lx/2, Lx/2))
quad1 = plt.pcolormesh(X, Y, Domain, cmap = cm.binary_r) 
scatter1 = plt.scatter(pos_x, pos_y, marker='o', s=1)

#Electric potential
ax2 = plt.subplot(gs[2, 0])
ax2.set(xlim=(-Lx/2, Lx/2), ylim=(-Lx/2, Lx/2))
quad2 = plt.pcolormesh(X, Y, Phi, cmap = cm.Blues) 
plt.colorbar()

#Electric field
ax3 = plt.subplot(gs[2, 1])
ax3.set(xlim=(-Lx/2, Lx/2), ylim=(-Lx/2, Lx/2))
stream1 = ax3.streamplot(U, V, Electric_field_x, Electric_field_y, linewidth=0.5, arrowsize=0.5, color='blue')
quad3 = plt.pcolormesh(X, Y, Electric_field, cmap = cm.Reds)
plt.colorbar()


#ax3.collections = [] # clear lines streamplot
#ax3.patches = [patch for patch in ax3.patches if not isinstance(x, mpl.patches.FancyArrowPatch)] # clear arrowheads streamplot

def init():
    quad1.set_array([])
    scatter1.set_offsets([])
    quad2.set_array([])
    stream1
    quad3.set_array([])
    return quad1,scatter1,quad2,stream1,quad3

def animate(Iter):
    #Load data
    n = Iter
    fig.suptitle('t = '+str(n*dt)+' s')
    fields = np.loadtxt('results/data/fields_'+str(n)+'.dat')
    particles = np.loadtxt('results/data/particles_'+str(n)+'.dat')

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

    #remove streamplot arrows
    for art in ax3.get_children():
        if not isinstance(art, mpl.patches.FancyArrowPatch):
            continue
        art.remove()        # Method 1


    quad1.set_array(Domain.ravel())
    scatter1.set_offsets(particles[:,:2])
    quad2.set_array(Phi.ravel())
    stream1 = ax3.streamplot(U, V, Electric_field_x, Electric_field_y, linewidth=0.5, arrowsize=0.5, color='blue')
    quad3.set_array(Electric_field.ravel())

    return quad1,scatter1,quad2,stream1.lines.remove(),quad3

gs.tight_layout(fig)
anim = animation.FuncAnimation(fig,animate,frames=frn,interval=interval,blit=False,repeat=False)
anim.save('animation.gif', writer='imagemagick', fps=fps)

#plt.show()
