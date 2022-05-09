import numpy as np
import matplotlib.pyplot as plt
import vtk


data=np.loadtxt('phi.dat')
x_index = data[:,0].astype('int32')
y_index = data[:,1].astype('int32')
z = data[:,2]

data=np.loadtxt('boundary.dat')
zz = data[:,2]

plt.imshow(z.reshape(max(x_index)+1,max(y_index)+1))
plt.show()

plt.imshow(zz.reshape(max(x_index)+1,max(y_index)+1))
plt.show()

n = max(x_index)+1
m = max(y_index)+1
xmin = 0; xmax = 10
ymin = 0; ymax = 10

plane = vtk.vtkPlaneSource()
plane.SetResolution(n-1,m-1)
plane.SetOrigin([xmin,ymin,0])  # Lower left corner
plane.SetPoint1([xmax,ymin,0])
plane.SetPoint2([xmin,ymax,0])
plane.Update()

# Map the values to the planar mesh.
# Assumption: same index i for scalars z[i] and mesh points
nPoints = plane.GetOutput().GetNumberOfPoints()
assert(nPoints == len(z))
# VTK has its own array format. Convert the input
# array (z) to a vtkFloatArray.
scalars = vtk.vtkFloatArray()
scalars.SetNumberOfValues(nPoints)
for i in range(nPoints):
    scalars.SetValue(i, z[i])
# Assign the scalar array.
plane.GetOutput().GetPointData().SetScalars(scalars)


writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName('output.vtp')
writer.SetInputConnection(plane.GetOutputPort())
writer.Write() # => Use for example ParaView to see scalars