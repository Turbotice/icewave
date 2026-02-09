# Code à executer dans un terminal pour avoir figure outline
import pyvista as pv
import numpy as np

# Load the SfePy result VTK file
filename = "C:/Users/Vasco Zanchi/Documents/git_local/sfepy_examples/its2D.vtk"
mesh = pv.read(filename)

# Check available point arrays
print("Point arrays in the mesh:", mesh.point_data.keys())

# Assuming displacement field is stored as 'u'
if 'u' not in mesh.point_data:
    raise ValueError("No displacement field 'u' found in the VTK file.")

displacements = mesh.point_data['u']

# Optionally, scale displacements for better visualization
scale_factor = 1.0  # adjust if deformation is too small/large
deformed_points = mesh.points + displacements * scale_factor

# Create a PyVista plotter
plotter = pv.Plotter()
plotter.add_mesh(mesh.points.copy(), show_edges=True, color='lightgrey')  # original mesh
plotter.add_mesh(deformed_points, color='orange', show_edges=True)          # deformed mesh
plotter.add_arrows(mesh.points, displacements, mag=1.0)                     # arrows for u
plotter.add_axes()
plotter.show_grid()
plotter.show()

# %%
