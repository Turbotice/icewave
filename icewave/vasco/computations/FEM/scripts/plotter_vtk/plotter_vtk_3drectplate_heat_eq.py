# %%
import pyvista as pv
import numpy as np
import glob
import matplotlib.pyplot as plt
from vtk.util.numpy_support import vtk_to_numpy

# -----------------------------
# PARAMÈTRES
# -----------------------------
disk = "D:"
results_dir = f"{disk}/FEM/results/time_dep_heat_eq"
pattern = results_dir + "/*.vtk"   # tous les fichiers VTK

# point dont on veut suivre la température
target_point = np.array([0.04, 0.02, 0.002])  # exemple (x,y,z) en m
tol = 1e-6  # tolérance pour comparaison

# -----------------------------
# 1. Lister et trier les fichiers VTK
# -----------------------------
filenames = sorted(glob.glob(pattern))
if not filenames:
    raise FileNotFoundError(f"Aucun fichier .vtk trouvé dans {results_dir}")

print(f"Trouvé {len(filenames)} fichiers VTK.")

# -----------------------------
# 2. Extraire la température au point choisi
# -----------------------------
times = []
temps = []

for i, fname in enumerate(filenames):
    mesh = pv.read(fname)

    # coordonnées du maillage
    points = vtk_to_numpy(mesh.GetPoints().GetData())

    # température scalaire
    if "T" not in mesh.point_data:
        raise ValueError(f"Pas de champ 'T' dans {fname}")

    Tfield = np.array(mesh.point_data["T"])

    # trouver le nœud le plus proche du point demandé
    dists = np.linalg.norm(points - target_point, axis=1)
    idx = np.argmin(dists)
    Tval = Tfield[idx]

    # estimer le temps depuis l'index (si pas stocké dans VTK)
    t = i  # tu peux adapter (par ex. i*dt si dt=0.5s)
    times.append(t)
    temps.append(Tval)

    print(f"{fname} -> T({target_point}) ≈ {Tval:.2f} K")

# -----------------------------
# 3. Plot T(t)
# -----------------------------
plt.figure()
plt.plot(times, temps, "o-", label=f"T @ {target_point}")
plt.xlabel("Time step")
plt.ylabel("Temperature [K]")
plt.title("Évolution de la température au point choisi")
plt.grid()
plt.legend()
plt.show()

