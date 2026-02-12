#%%
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy

# -----------------------------
# 1. Load the SfePy result
# -----------------------------
disk = 'D:'


filename = f"{disk}/FEM/results/plate_E3.0GPa_L80.0mm_w40.0mm_h4.0mm_Fz-10.0N.vtk"


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


# Exemple d'utilisation :
filename = make_filename(E_GPa=3.0, L_mm=80.0, w_mm=40.0, h_mm=4.0, Fz=-30)
filepath = f"{disk}/FEM/results/{filename}"
print(filepath)

mesh = pv.read(filepath)

# -----------------------------
# 2. Extract points and displacements
# -----------------------------
points = vtk_to_numpy(mesh.GetPoints().GetData())  # Nx3
displacements = np.array(mesh.point_data['u'])

# Pad 2D → 3D if needed
if displacements.shape[1] == 2:
    displacements = np.hstack([displacements, np.zeros((displacements.shape[0], 1))])

xyz = points[:, :3]
u = displacements[:, :3]

# -----------------------------
# 3. Select line (y ≈ y0, z ≈ z0)
# -----------------------------
y0 = np.mean(xyz[:,1])   # milieu en y
z0 = np.min(xyz[:,2])    # couche inférieure
tol = 1e-6

mask = (np.abs(xyz[:,1] - y0) < tol) & (np.abs(xyz[:,2] - z0) < tol)

x_sel = xyz[mask, 0]
uz_sel = u[mask, 2]

# -----------------------------
# 4. Plot profile u_z(x)
# -----------------------------
order = np.argsort(x_sel)
plt.figure()
plt.plot(x_sel[order], uz_sel[order], 'o-')
plt.xlabel("x [m]")
plt.ylabel("u_z [m]")
plt.title("Deflection along bottom middle line")
plt.grid()
plt.show()

# %%
