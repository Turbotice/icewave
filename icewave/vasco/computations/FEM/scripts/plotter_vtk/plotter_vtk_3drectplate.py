#%%
import pyvista as pv
import numpy as np
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib.pyplot as plt

# -----------------------------
# 1. Load the SfePy result
# -----------------------------
disk = 'E:'

def make_filename(E_GPa, L_mm, w_mm, h_mm, Fz):
    """
    Génère un nom de fichier de type :
    plate_E3.0GPa_L80.0mm_w40.0mm_h4.0mm_Fz-10.0N.vtk
    """
    return (f"plate_E{E_GPa:.1f}GPa"
            f"_L{L_mm:.0f}mm"
            f"_w{w_mm:.0f}mm"
            f"_h{h_mm:.2f}mm"
            f"_Fz{Fz:.1f}N.vtk")

filename = make_filename(3, 80, 40,4, -10) #f"{disk}/FEM/results/plate_E3.0GPa_L80.0mm_w40.0mm_h4.0mm_Fz-10.0N.vtk"
filepath = f"{disk}/FEM/results/{filename}"
mesh = pv.read(filepath)

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
scale_factor = 1
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

# c) y component displacement
vertical_disp = displacements[:, 1]
pv_mesh["u_y"] = vertical_disp

# c) vertical displacement (z component)
if displacements.shape[1] == 3:
    vertical_disp_z = displacements[:, 2]  # u_z
else:
    vertical_disp_z = np.zeros(displacements.shape[0])  # si 2D

pv_mesh["u_z"] = vertical_disp_z


# -----------------------------
# 6. Plot
# -----------------------------
plotter = pv.Plotter(shape=(1, 1))  

# Left: displacement magnitude
plotter.subplot(0, 0)
plotter.add_text("Displacement magnitude", font_size=12)
plotter.add_mesh(pv_mesh, scalars="disp_magnitude", show_edges=True, cmap="jet")
plotter.add_axes()
plotter.show_grid()
plotter.show()

# Bottom left: horizontal displacement
plotter = pv.Plotter(shape=(1, 1))  

plotter.subplot(0, 0)
plotter.add_text("Horizontal displacement (u_x)", font_size=12)
plotter.add_mesh(pv_mesh, scalars="u_x", show_edges=True, cmap="viridis")
plotter.add_axes()
plotter.show_grid()
plotter.show()

# Bottom Right: vertical displacement
plotter = pv.Plotter(shape=(1, 1))  

plotter.subplot(0, 0)
plotter.add_text("Vertical displacement (u_y)", font_size=12)
plotter.add_mesh(pv_mesh, scalars="u_y", show_edges=True, cmap="coolwarm")
plotter.add_axes()
plotter.show_grid()
#plotter.link_views()  # link camera between subplots (si on veut tout voir en même temps)
plotter.show()

plotter = pv.Plotter(shape=(1,1))
plotter.add_mesh(pv_mesh, scalars="u_z", show_edges=True, cmap="coolwarm")
plotter.add_text("Vertical displacement u_z")
plotter.add_axes()
plotter.show_grid()
plotter.show()

# %% computation of epsilon_xx and epsilon_xy


import numpy as np
from scipy.spatial import KDTree

xyz = points[:, :3]
u = displacements[:, :3]

tree = KDTree(xyz)
k = 10  # number of neighbors

# Initialize strain arrays
eps_xx = np.zeros(len(xyz))
eps_yy = np.zeros(len(xyz))
eps_zz = np.zeros(len(xyz))
eps_xy = np.zeros(len(xyz))
eps_xz = np.zeros(len(xyz))
eps_yz = np.zeros(len(xyz))

for i, (xi, yi, zi) in enumerate(xyz):
    dists, idxs = tree.query([xi, yi, zi], k=k)
    neighbors = xyz[idxs]
    u_neighbors = u[idxs]

    # Fit planes u_x = a*x + b*y + c*z + d
    A = np.c_[neighbors, np.ones(len(neighbors))]
    coeff_x, _, _, _ = np.linalg.lstsq(A, u_neighbors[:,0], rcond=None)
    coeff_y, _, _, _ = np.linalg.lstsq(A, u_neighbors[:,1], rcond=None)
    coeff_z, _, _, _ = np.linalg.lstsq(A, u_neighbors[:,2], rcond=None)

    # Small strains
    eps_xx[i] = coeff_x[0]
    eps_yy[i] = coeff_y[1]
    eps_zz[i] = coeff_z[2]
    eps_xy[i] = 0.5*(coeff_x[1] + coeff_y[0])
    eps_xz[i] = 0.5*(coeff_x[2] + coeff_z[0])
    eps_yz[i] = 0.5*(coeff_y[2] + coeff_z[1])


pv_mesh["eps_xx"] = eps_xx
pv_mesh["eps_yy"] = eps_yy
pv_mesh["eps_zz"] = eps_zz
pv_mesh["eps_xy"] = eps_xy
pv_mesh["eps_xz"] = eps_xz
pv_mesh["eps_yz"] = eps_yz

plotter = pv.Plotter(shape=(2,3))

plotter.subplot(0,0)
plotter.add_mesh(pv_mesh, scalars="eps_xx", cmap="jet", show_edges=True)
plotter.add_text("ε_xx")

plotter.subplot(0,1)
plotter.add_mesh(pv_mesh, scalars="eps_yy", cmap="jet", show_edges=True)
plotter.add_text("ε_yy")

plotter.subplot(0,2)
plotter.add_mesh(pv_mesh, scalars="eps_zz", cmap="jet", show_edges=True)
plotter.add_text("ε_zz")

plotter.subplot(1,0)
plotter.add_mesh(pv_mesh, scalars="eps_xy", cmap="coolwarm", show_edges=True)
plotter.add_text("ε_xy")

plotter.subplot(1,1)
plotter.add_mesh(pv_mesh, scalars="eps_xz", cmap="coolwarm", show_edges=True)
plotter.add_text("ε_xz")

plotter.subplot(1,2)
plotter.add_mesh(pv_mesh, scalars="eps_yz", cmap="coolwarm", show_edges=True)
plotter.add_text("ε_yz")

plotter.link_views()
plotter.show()


# %%
y0 = np.mean(xyz[:,1])  # milieu en y
z0 = np.min(xyz[:,2])   # surface supérieure
tol = 1e-6

mask = (np.abs(xyz[:,1] - y0) < tol) & (np.abs(xyz[:,2] - z0) < tol)

x_sel = xyz[mask, 0]
uz_sel = u[mask, 2]

order = np.argsort(x_sel)
plt.figure()
plt.plot(x_sel[order], uz_sel[order], 'o-')
plt.xlabel("x [m]")
plt.ylabel("u_z [m]")
plt.title("Deflection along top middle line")
plt.grid()
plt.show()

# %%
