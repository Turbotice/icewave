#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%% create synthetic Finite Elements data :
# Simple mesh
x = np.linspace(-5, 5, 21)
y = np.linspace(-5, 5, 21)
X, Y = np.meshgrid(x, y)

# Fake displacements (mode I style)
ux = 0.001 * X
uy = -0.001 * Y

# Build nodemap structure (facet_id, x, y, z, ux, uy, uz)
df = pd.DataFrame({
    "facet_id": np.arange(X.size),
    "x": X.flatten(),
    "y": Y.flatten(),
    "z": 0.0,             # flat 2D geometry
    "ux": ux.flatten(),
    "uy": uy.flatten(),
    "uz": 0.0,            # no out-of-plane displacement
})

FE_folder = "C:/Users/Vasco Zanchi/Documents/crack_analysis/FE_data/"
fname = "synthetic_fe_nodemap.txt"

# Save as a text file (mimics nodemap)
df.to_csv(FE_folder + fname, sep="\t", index=False)



# %%
"""
from crackpy.fracture_analysis.data_processing import InputData
from crackpy.structure_elements.data_files import Nodemap
from crackpy.fracture_analysis.line_integration import *

# 1. Prepare material
material = Material(E=210000, nu_xy=0.3, plane_strain=True)

# 2. Load a nodemap (replace with your real path + filename)
nodemap = Nodemap(folder=FE_folder, name="synthetic_fe_nodemap.txt")

# 3. Wrap into InputData
data = InputData(nodemap)
data.calc_stresses()   # compute stresses from displacements

# 4. Define integration path
path_props = PathProperties(size_left=-1, size_right=1,
                            size_bottom=-1, size_top=1,
                            tick_size=0.1,
                            top_offset=0.2, bottom_offset=-0.2)
int_path = IntegrationPath(origin_x=0.0, origin_y=0.0, path_properties=path_props)
int_path.create_nodes()

# 5. Run LineIntegral
line_int = LineIntegral(int_path, data, material)
line_int.integrate()

# 6. Results
print("J =", line_int.j_integral)
print("K_J =", line_int.sif_k_j)
"""
#%%
import numpy as np
import pandas as pd

# Maillage
x = np.linspace(-5, 5, 21)
y = np.linspace(-5, 5, 21)
X, Y = np.meshgrid(x, y)

# Crack tip au centre
x_tip, y_tip = 0.0, 0.0

# Distance au crack tip
r = np.sqrt((X - x_tip)**2 + (Y - y_tip)**2)
r[r==0] = 1e-6  # éviter division par zéro

# Amplitude du déplacement
A = 0.01

# Déplacements Mode I synthétiques
ux = 0.3 * A * np.sqrt(r)      # déplacement horizontal faible
uy = A * np.sqrt(r)            # déplacement vertical principal
uz = np.zeros_like(ux)         # pas de déplacement hors plan

# Construire DataFrame type FE nodemap
df = pd.DataFrame({
    "facet_id": np.arange(X.size),
    "x": X.flatten(),
    "y": Y.flatten(),
    "z": 0.0,
    "ux": ux.flatten(),
    "uy": uy.flatten(),
    "uz": uz.flatten(),
})

# Sauvegarde
FE_folder = "C:/Users/Vasco Zanchi/Documents/crack_analysis/FE_data/"
fname = "synthetic_fe_modeI.txt"
df.to_csv(FE_folder + fname, sep="\t", index=False)

print(df.head())



#%%




import numpy as np
from crackpy.fracture_analysis.data_processing import InputData
from crackpy.structure_elements.material import Material

# Load your tab-separated FE nodemap
raw = np.genfromtxt(FE_folder+fname, delimiter="\t", skip_header=1)

facet_id = raw[:, 0]
coor_x = raw[:, 1]
coor_y = raw[:, 2]
coor_z = raw[:, 3]
disp_x = raw[:, 4]
disp_y = raw[:, 5]
disp_z = raw[:, 6]

# If you don't have strains, you can compute approximate strains from displacements, or set dummy values
eps_x = np.zeros_like(coor_x)
eps_y = np.zeros_like(coor_x)
eps_xy = np.zeros_like(coor_x)

# Initialize InputData manually
data = InputData()
data.set_data_manually(coor_x, coor_y, disp_x, disp_y, eps_x, eps_y, eps_xy)

# Define material
material = Material(E=210000, nu_xy=0.3, plane_strain=True)

# Compute stresses from displacements
data.calc_stresses(material)
data.calc_eps_vm()

print(data.eps_vm[:5])
print(data.sig_vm[:5])



# %%

from crackpy.fracture_analysis.line_integration import PathProperties, IntegrationPath, LineIntegral

# Crack tip au centre
x_tip, y_tip = 0.0, 0.0

# Définir les propriétés du chemin
path_props = PathProperties(
    size_left=-1, size_right=1,
    size_bottom=-1, size_top=1,
    tick_size=0.1
)

# Créer le chemin
int_path = IntegrationPath(origin_x=x_tip, origin_y=y_tip, path_properties=path_props)
int_path.create_nodes()

# %%
