#%%
import pyvista as pv
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy

# -----------------------------
# 1. Load the SfePy result
# -----------------------------
filename = "C:/Users/Vasco Zanchi/Documents/git_local/sfepy_examples/brazilian_test/its2D.vtk"
mesh = pv.read(filename)

# -----------------------------
# 2. Convert points and displacement vectors to NumPy
# -----------------------------
points = vtk_to_numpy(mesh.GetPoints().GetData())  # Nx3 array

# Displacements
if 'u' not in mesh.point_data:
    raise ValueError("No displacement field 'u' found in the VTK file.")
displacements = np.array(mesh.point_data['u'])   # Nx2 or Nx3

# Pad 2D displacements to 3D
if displacements.shape[1] == 2:
    displacements = np.hstack([displacements, np.zeros((displacements.shape[0], 1))])

# -----------------------------
# 3. Compute deformed points
# -----------------------------
scale_factor = 1.0
deformed_points = points + displacements * scale_factor

# -----------------------------
# 4. Use the original PyVista mesh, update points
# -----------------------------
pv_mesh = mesh.copy()
pv_mesh.points = deformed_points

# -----------------------------
# 5. Add colormaps
# -----------------------------
# a) displacement magnitude (for reference)
disp_magnitude = np.linalg.norm(displacements, axis=1)
pv_mesh["disp_magnitude"] = disp_magnitude

# b) horizontal displacement (x component)
horizontal_disp = displacements[:, 0]
pv_mesh["u_x"] = horizontal_disp

# c) vertical displacement (y component)
vertical_disp = displacements[:, 1]
pv_mesh["u_y"] = vertical_disp

# -----------------------------
# 6. Plot
# -----------------------------
plotter = pv.Plotter(shape=(1, 1))  # side-by-side plots

# Left: displacement magnitude
plotter.subplot(0, 0)
plotter.add_text("Displacement magnitude", font_size=12)
plotter.add_mesh(pv_mesh, scalars="disp_magnitude", show_edges=True, cmap="jet")
plotter.add_axes()
plotter.show_grid()
plotter.show()

# Bottom left: horizontal displacement
plotter = pv.Plotter(shape=(1, 1))  # side-by-side plots

plotter.subplot(0, 0)
plotter.add_text("Horizontal displacement (u_x)", font_size=12)
plotter.add_mesh(pv_mesh, scalars="u_x", show_edges=True, cmap="viridis")
plotter.add_axes()
plotter.show_grid()
plotter.show()

# Bottom Right: vertical displacement
plotter = pv.Plotter(shape=(1, 1))  # side-by-side plots

plotter.subplot(0, 0)
plotter.add_text("Vertical displacement (u_y)", font_size=12)
plotter.add_mesh(pv_mesh, scalars="u_y", show_edges=True, cmap="coolwarm")
plotter.add_axes()
plotter.show_grid()

#plotter.link_views()  # link camera between subplots (si on veut tout voir en même temps)

plotter.show()


# %% computation of epsilon_xx and epsilon_xy

import numpy as np
from scipy.spatial import KDTree

# 2D points and displacements
xy = points[:, :2]           # Nx2
u = displacements[:, :2]     # Nx2

# Build KDTree for nearest neighbors
tree = KDTree(xy)

# Initialize strain arrays
eps_xx = np.zeros(xy.shape[0])
eps_xy = np.zeros(xy.shape[0])

# Loop over each node
for i, (xi, yi) in enumerate(xy):
    # find nearest neighbors (e.g., 4 closest nodes)
    dists, idxs = tree.query([xi, yi], k=5)
    neighbors = xy[idxs]
    u_neighbors = u[idxs]

    # Fit a plane to u_x = a*x + b*y + c
    A = np.c_[neighbors, np.ones(len(neighbors))]
    coeff_x, _, _, _ = np.linalg.lstsq(A, u_neighbors[:,0], rcond=None)
    coeff_y, _, _, _ = np.linalg.lstsq(A, u_neighbors[:,1], rcond=None)

    # Strain components
    eps_xx[i] = coeff_x[0]               # du_x/dx
    eps_xy[i] = 0.5*(coeff_x[1] + coeff_y[0])  # 0.5*(du_x/dy + du_y/dx)



pv_mesh["eps_xx"] = eps_xx
pv_mesh["eps_xy"] = eps_xy

plotter = pv.Plotter(shape=(1,2))

plotter.subplot(0,0)
plotter.add_mesh(pv_mesh, scalars="eps_xx", cmap="jet", show_edges=True)
plotter.add_text("Strain ε_xx")

plotter.subplot(0,1)
plotter.add_mesh(pv_mesh, scalars="eps_xy", cmap="coolwarm", show_edges=True)
plotter.add_text("Strain ε_xy")

plotter.link_views()
plotter.show()

# %%
## Pour etre sûr faire un cercle pour le test brésilien...  