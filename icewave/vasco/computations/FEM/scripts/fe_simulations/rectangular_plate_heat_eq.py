#%%
"""
Transient heat conduction in a plate:
- T prescribed on bottom face
- Heat flux prescribed on top face
"""

from __future__ import absolute_import
import numpy as np
from sfepy import data_dir
from sfepy.base.base import Struct

disk = 'D:'
data_dir = f'{disk}/FEM/'

# -----------------------------
# INPUTS
# -----------------------------
T_bottom = float(input("Température imposée sur la face inférieure (°C) ? "))
q_top = float(input("Flux imposé sur la face supérieure (W/m²) ? "))
rho = float(input("Masse volumique (kg/m³) ? "))
cp = float(input("Capacité thermique (J/kg/K) ? "))
k_val = float(input("Conductivité thermique (W/m/K) ? "))

filename_mesh = data_dir + '/meshes/3d/3drect_plate.mesh'
filename_info = data_dir + '/meshes/3d/plate_info.txt'

plate_info = np.loadtxt(filename_info)
plate_length_x = plate_info[0]
plate_width_y = plate_info[1]
plate_thickness_z = plate_info[2]

tol = 1e-6

# -----------------------------
# Post-processing
# -----------------------------
def post_proces(out, pb, state, extend=False):
    # Sauvegarde avec le temps courant dans le nom du fichier
    filename = (f"{pb.output_dir}/heatPlate_"
                f"Tbot{T_bottom:.1f}C_qtop{q_top:.1f}_"
                f"t{pb.ts.time:.3f}.vtk")
    pb.save_state(filename, state, out=out)
    print(f"Saved: {filename}")
    return out

# -----------------------------
# Options
# -----------------------------
options = {
    'post_process_hook': 'post_proces',
    'output_dir': f'{disk}/FEM/results/time_dep_heat_eq',
    'output_format': None,  # désactive sauvegarde auto, on gère dans post_proces
}

# -----------------------------
# Regions
# -----------------------------
regions = {
    'Omega': 'all',
    'Bottom': (f'vertices in (z < {0.0 + tol})', 'facet'),
    'Top':    (f'vertices in (z > {plate_thickness_z - tol})', 'facet'),
}

# -----------------------------
# Materials
# -----------------------------
materials = {
    'plate': ({
        'lam': k_val,
        'rho_cp': rho * cp,   # effective volumetric heat capacity
    },),
    'flux_top': ({'val': q_top},),
}


# -----------------------------
# Fields
# -----------------------------
fields = {
    'temperature': ('real', 1, 'Omega', 1),
}

# -----------------------------
# Integrals
# -----------------------------
integrals = {
    'i': 2,
}

# -----------------------------
# Variables
# -----------------------------
variables = {
    'T': ('unknown field', 'temperature', 0, 1),
    's': ('test field', 'temperature', 'T'),
}

# -----------------------------
# Boundary conditions
# -----------------------------
ebcs = {
    'T_bottom': ('Bottom', {'T.0': T_bottom}),
}

# -----------------------------
# Equations
# -----------------------------
equations = {
    'heat_balance':
    """
    dw_dot.i.Omega(plate.rho_cp, s, T)
    + dw_laplace.i.Omega(plate.lam, s, T)
    = dw_surface_integrate.i.Top(flux_top.val, s)
    """
}


# -----------------------------
# Solvers
# -----------------------------
solvers = {
    'ls': ('ls.auto_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
        'eps_a': 1e-8,
    }),
    'ts': ('ts.simple', {
        't0': -10,    # temps initial
        't1': 10.0,   # temps final (s)
        'dt': 0.1,    # pas de temps
        'n_step': None,  # calculé automatiquement via t0, t1, dt
    }),
}
