#%%
r"""
A linear elastic beam loaded with a continuous force. The FE mesh consists
of hexehedral and tetrahedral elements.

The displacement at the beam end is compared to the reference
solution calculated on the homogeneous hexahedral mesh.

Running the simulation::

    sfepy-run sfepy/examples/linear_elasticity/mixed_mesh.py

Viewing the results::

    sfepy-view beam_h* -f u:s0:wu:e:p0 u:s1:wu:e:p0 --camera-position="1.2,-0.6,0.1,0.4,0.1,-0.1,-0.2,0.1,1"
"""
from __future__ import absolute_import
from sfepy import data_dir
from sfepy.mechanics.matcoefs import stiffness_from_youngpoisson
from sfepy.base.base import Struct
import numpy as np
from datetime import datetime

disk = 'D:'
data_dir = f'{disk}/FEM/'


def get_force(ts, coors, mode=None, **kwargs): # fonction pas utilisée
    if mode == 'qp':
        F = 1e3

        val = np.zeros_like(coors)[..., None]
        val[:, 2, 0] = -coors[:, 0] / 0.7 * F

        return {'val': val}



'''
def post_proces(out, pb, state, extend=False):
    c1 = pb.domain.regions['Omega'].cells
    #c2 = pb.domain.regions['Omega_t'].cells
    ed = pb.domain.regions['Edge'].vertices

    S = np.zeros((pb.domain.cmesh.n_el, 1, 6, 1), dtype=np.float64)

    S[c1] = pb.evaluate('ev_cauchy_stress.i.Omega(solid.D, u)',
                        mode='el_avg')
    #S[c2] = pb.evaluate('ev_cauchy_stress.i.Omega_t(solid.D, u_t)',
    #                    mode='el_avg')

    out['stress'] = Struct(name='out', data=S, mode='cell')

    u_nd = out['u'].data[ed[0], :]

    #ref_u_nd = reference_solution(pb).reshape(-1, 3)[ed[0], :]
    #du = np.linalg.norm((ref_u_nd - u_nd) / np.linalg.norm(ref_u_nd))
    #print(f'Relative difference with respect to the reference solution: {du}')

    return out

'''


def post_proces(out, pb, state, extend=False):
    c1 = pb.domain.regions['Omega'].cells
    ed = pb.domain.regions['Edge'].vertices

    S = np.zeros((pb.domain.cmesh.n_el, 1, 6, 1), dtype=np.float64)
    S[c1] = pb.evaluate('ev_cauchy_stress.i.Omega(solid.D, u)', mode='el_avg')

    out['stress'] = Struct(name='out', data=S, mode='cell')

    u_nd = out['u'].data[ed[0], :]

    # ---- custom filename ----
    #filename = f"{pb.output_dir}/plate_E{np.round(E/1e9,2)}GPa_L{np.round(plate_length_x*1e3,0)}mm_w{np.round(plate_width_y*1e3,0)}mm_h{np.round(plate_thickness_z*1e3,2)}mm_Fz{np.round(Fz,2)}N.vtk"
    filename = (
        f"{pb.output_dir}/plate_"
        f"E{np.round(E/1e9,2)}GPa_"          # 2 digits après la virgule
        f"L{plate_length_x*1e3:.0f}mm_"  # entier (pas de décimales)
        f"w{plate_width_y*1e3:.0f}mm_"   # entier
        f"h{plate_thickness_z*1e3:.2f}mm_"  # 2 digits
        f"Fz{Fz:.1f}N.vtk"           # 2 digits
        )

    pb.save_state(filename, state, out=out)
    print(f"Saved custom file: {filename}")

    return out



##############################################################
# INPUTS
#Fz = -10 # force along z in Newtons (for 3 points bending, negative sign)
Fz = -0.5 # float(input('Applied force Fz (N) ? '))
E = 0.03e9 # float(input("Young's modulus (GPa) ? ")) * 1e9

plate_length_x = 8e-2
plate_width_y = 4e-2
plate_thickness_z = 9.7e-3
nx_cells = 40
ny_cells = 20
nz_cells = 10

filename_mesh = f'{data_dir}/meshes/3d/3drect_plate_L{np.round(plate_length_x*1e3).astype(int)}mm_w{np.round(plate_width_y*1e3).astype(int)}mm_h{np.round(plate_thickness_z*1e3).astype(int)}mm_nx{nx_cells}_ny{ny_cells}_nz{nz_cells}.mesh'
filename_nodesTopMiddleLine = f'{data_dir}/meshes/3d/nodes_TopMiddleLine_L{np.round(plate_length_x*1e3).astype(int)}mm_w{np.round(plate_width_y*1e3).astype(int)}mm_h{np.round(plate_thickness_z*1e3).astype(int)}mm_nx{nx_cells}_ny{ny_cells}_nz{nz_cells}.txt'


data_nodesTopMiddleLine = np.loadtxt(filename_nodesTopMiddleLine)
arr_indicesTopMiddleLine = data_nodesTopMiddleLine[:,0]
str_indices_TopMiddleLine = ','.join(str(int(x)) for x in arr_indicesTopMiddleLine)
###############################################################

options = {
    'post_process_hook': 'post_proces',
    'output_dir': f'{disk}/FEM/results/',  # <--- folder where vtk files will be saved
    'output_format': None,  # disable automatic saving
}


tol = 1e-6  # or slightly larger if needed

regions = {

    # whole domain
    'Omega_': 'all',

    # groups for hexas and tetras (if you still output both)
    'Omega':   ('cells of group 1', 'cell', None, {'vertices_from': 'Omega_'}),
    #'Omega_t': ('cells of group 2', 'cell', None, {'vertices_from': 'Omega_'}),

    # supports at left edge (x ~ 0)
    'Left':  (f'vertices in (x < {0.0 + tol})', 'facet'),
    # supports at right edge (x ~ Lx)
    'Right': (f'vertices in (x > {plate_length_x - tol})', 'facet'),

    # top and bottom surfaces in z
    'Top':   (f'vertices in (z > {plate_thickness_z - tol})', 'facet'),
    'Bottom':(f'vertices in (z < {0.0 + tol})', 'facet'),

    # bottom left
    'BottomLeft':(
        f'vertices in (z < {0.0 + tol}) & '
        f'(x < {0.0 + tol})',
         'vertex'),
    # bottom left
    'BottomRight':(
        f'vertices in (z < {0.0 + tol}) & '
        f'(x > {plate_length_x - tol})',
         'vertex'),

    'Edge': (
        f'vertices in (x > {plate_length_x - tol}) & (z < {0.0 + tol})',
        'vertex'),
    'TopMiddleLine': (
        f'vertex {str_indices_TopMiddleLine}','vertex'
    ),

}

"""
    # optional: edges for loading/measurement
    'TopMiddleLine': (
        f'vertices in (x > {0.49*plate_length_x - tol}) & '
        f'(x < {0.51*plate_length_x + tol}) & '
        f'(z > {plate_thickness_z - tol})',
        'vertex'
    ),
"""
functions = {
    'get_force' : (get_force,),
}



materials = {
    'solid': ({'D': stiffness_from_youngpoisson(dim=3, young=E, poisson=0.3)},),
    'force': 'get_force',
    'vforce': ({'.val' : [0.0, 0.0, Fz/len(arr_indicesTopMiddleLine)]},),
}

fields = {
    'displacement': ('real', 'vector', 'Omega', 1),
    #'displacement_t': ('real', 'vector', 'Omega_t', 1),
}

integrals = {
    'i': 2,
}

variables = {
    'u': ('unknown field', 'displacement', 0),
    'v': ('test field', 'displacement', 'u'),
    #'u_t': ('unknown field', 'displacement_t', 0),
    #'v_t': ('test field', 'displacement_t', 'u_t'),
}

"""
ebcs = {
    'FixedBottomLeft': ('BottomLeft', {'u.all' : 0.0}),
    'FixedBottomRight': ('BottomRight', {'u.all' : 0.0}),
}
"""
# correction des conditions limites : on impose seulement la composante verticale de la vitesse nulle dans les coins inférieurs
ebcs = {
    'FixedBottomLeft': ('BottomLeft', {'u.2' : 0.0}),
    'FixedBottomRight': ('BottomRight', {'u.2' : 0.0}),
}


equations = {
    'balance_of_forces':
    """dw_lin_elastic.i.Omega(solid.D, v, u)
     = dw_point_load.0.TopMiddleLine(vforce.val, v)
     """
}



solvers = {
    'ls': ('ls.auto_direct', {}),
    'newton': ('nls.newton', {
        'i_max': 1,
        'eps_a': 1e-6,
    }),
}


# %%
