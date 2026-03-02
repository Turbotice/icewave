#%%
import numpy as np
import pyvista as pv
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

def read_avg_plate_profile(E_GPa=3.0, L_mm=80.0, w_mm=40.0, h_mm=4.0, Fz=-30, plot=False):
    # Exemple d'utilisation :
    filename = make_filename(E_GPa=E_GPa, L_mm=L_mm, w_mm=w_mm, h_mm=h_mm, Fz=Fz)
    filepath = f"{disk}/FEM/results/{filename}"
    print(filepath)

    # Chargement du maillage
    mesh = pv.read(filepath)

    # Extraction
    points = vtk_to_numpy(mesh.GetPoints().GetData())  # Nx3
    displacements = np.array(mesh.point_data["u"])

    # Ajustement 2D -> 3D si nécessaire
    if displacements.shape[1] == 2:
        displacements = np.hstack([displacements, np.zeros((displacements.shape[0], 1))])

    xyz = points[:, :3]
    u = displacements[:, :3]

    x = xyz[:, 0]
    uz = u[:, 2]

    # ------------------------------------------
    # 1. Regroupement des points selon x
    # ------------------------------------------

    # Tolérance : regrouper les x presque identiques
    tol = 1e-6

    # Arrondir les x pour créer des "bins"
    x_round = np.round(x / tol) * tol

    # Trouver les valeurs uniques
    x_unique = np.unique(x_round)

    x_mean = []
    uz_mean = []
    uz_std = []

    for xu in x_unique:
        mask = (x_round == xu)
        x_mean.append(np.mean(x[mask]))
        uz_mean.append(np.mean(uz[mask]))
        uz_std.append(np.std(uz[mask]))

    x_mean = np.array(x_mean)
    uz_mean = np.array(uz_mean)
    uz_std = np.array(uz_std)

    # Trier
    order = np.argsort(x_mean)

    # ------------------------------------------
    # 2. Plot avec écart-type
    # ------------------------------------------
    if plot:
        plt.figure()
        plt.plot(x_mean[order], uz_mean[order], '-', label='Déflexion moyenne')
        plt.fill_between(
            x_mean[order],
            uz_mean[order] - uz_std[order],
            uz_mean[order] + uz_std[order],
            alpha=0.3,
            label='Écart-type'
        )

        plt.xlabel("x [m]")
        plt.ylabel("u_z [m]")
        plt.title("Déflexion moyenne u_z(x) avec écart-type (moyenne en y,z)")
        plt.grid()
        plt.legend()
        plt.show()
    
    return x_mean[order], uz_mean[order], uz_std[order]

# %%
E_GPa = 3.0
L_mm = 80.0
w_mm = 40.0
h_mm=4.0

array_Fz = np.array([2,4,6,8,10,12,14,16,18,20],dtype=float) * (-1)
#array_Fz = np.array([2,4,6,8,10],dtype=float) * (-1)

array_fleche_mean = np.zeros(len(array_Fz))
array_fleche_std = np.zeros(len(array_Fz))

dict_results = {}
for i in range(len(array_Fz)):
    key_simu = make_filename(E_GPa, L_mm, w_mm, h_mm, array_Fz[i])
    dict_results[key_simu] = {}
    dict_results[key_simu]['xmean'], dict_results[key_simu]['uz_mean'], dict_results[key_simu]['uz_std'] = read_avg_plate_profile(E_GPa=E_GPa, L_mm=L_mm, w_mm=w_mm, h_mm=h_mm, Fz=array_Fz[i])
    argmax = np.argmax(np.abs(dict_results[key_simu]['uz_mean']))
    print(argmax)
    dict_results[key_simu]['uz_max_mean'] = dict_results[key_simu]['uz_mean'][argmax]
    dict_results[key_simu]['uz_max_std'] = dict_results[key_simu]['uz_std'][argmax]
    
    array_fleche_mean[i] = dict_results[key_simu]['uz_max_mean']
    array_fleche_std[i] = dict_results[key_simu]['uz_max_std']

# experssion théorique (poutres) de Fz vs deflection au centre
array_fleche_th = (1/(4*w_mm*1e-3*E_GPa*1e9)) * ((L_mm/h_mm)**3) * array_Fz


plt.figure()
plt.errorbar(np.abs(array_Fz), np.abs(array_fleche_mean), array_fleche_std, linestyle='', ecolor='k',color='b',marker='.',label='FEM simulations')
plt.plot(np.abs(array_Fz), np.abs(array_fleche_th),color='tab:orange',label='elastic beam theory')
#plt.xlim(0, np.max(np.abs(array_Fz))*1.1)
#plt.ylim(0, np.max(np.abs(array_fleche_mean))*1.1)
plt.xlabel('Fz magnitude (N)')
plt.ylabel('Deflection at the center (m)')
plt.legend()


# et pour 9GPa

E_GPa = 9.0
L_mm = 80.0
w_mm = 40.0
h_mm=4.0

array_Fz = np.array([2,4,6,8,10,12,14,16,18,20],dtype=float) * (-1)
#array_Fz = np.array([2,4,6,8,10],dtype=float) * (-1)

array_fleche_mean = np.zeros(len(array_Fz))
array_fleche_std = np.zeros(len(array_Fz))

for i in range(len(array_Fz)):
    key_simu = make_filename(E_GPa, L_mm, w_mm, h_mm, array_Fz[i])
    dict_results[key_simu] = {}
    dict_results[key_simu]['xmean'], dict_results[key_simu]['uz_mean'], dict_results[key_simu]['uz_std'] = read_avg_plate_profile(E_GPa=E_GPa, L_mm=L_mm, w_mm=w_mm, h_mm=h_mm, Fz=array_Fz[i])
    argmax = np.argmax(np.abs(dict_results[key_simu]['uz_mean']))
    print(argmax)
    dict_results[key_simu]['uz_max_mean'] = dict_results[key_simu]['uz_mean'][argmax]
    dict_results[key_simu]['uz_max_std'] = dict_results[key_simu]['uz_std'][argmax]
    
    array_fleche_mean[i] = dict_results[key_simu]['uz_max_mean']
    array_fleche_std[i] = dict_results[key_simu]['uz_max_std']


# experssion théorique (poutres) de Fz vs deflection au centre
array_fleche_th = (1/(4*w_mm*1e-3*E_GPa*1e9)) * ((L_mm/h_mm)**3) * array_Fz

plt.errorbar(np.abs(array_Fz), np.abs(array_fleche_mean), array_fleche_std, linestyle='', ecolor='k',color='b',marker='.',label='FEM simulations')
plt.plot(np.abs(array_Fz), np.abs(array_fleche_th),color='tab:orange',label='elastic beam theory')



# et pour epaisseur differente

E_GPa = 3.0
L_mm = 80.0
w_mm = 40.0
h_mm=2.0

array_Fz = np.array([2,4,6,8,10],dtype=float) * (-1)
#array_Fz = np.array([2,4,6,8,10],dtype=float) * (-1)

array_fleche_mean = np.zeros(len(array_Fz))
array_fleche_std = np.zeros(len(array_Fz))

for i in range(len(array_Fz)):
    key_simu = make_filename(E_GPa, L_mm, w_mm, h_mm, array_Fz[i])
    dict_results[key_simu] = {}
    dict_results[key_simu]['xmean'], dict_results[key_simu]['uz_mean'], dict_results[key_simu]['uz_std'] = read_avg_plate_profile(E_GPa=E_GPa, L_mm=L_mm, w_mm=w_mm, h_mm=h_mm, Fz=array_Fz[i])
    argmax = np.argmax(np.abs(dict_results[key_simu]['uz_mean']))
    print(argmax)
    dict_results[key_simu]['uz_max_mean'] = dict_results[key_simu]['uz_mean'][argmax]
    dict_results[key_simu]['uz_max_std'] = dict_results[key_simu]['uz_std'][argmax]
    
    array_fleche_mean[i] = dict_results[key_simu]['uz_max_mean']
    array_fleche_std[i] = dict_results[key_simu]['uz_max_std']

# experssion théorique (poutres) de Fz vs deflection au centre
array_fleche_th = (1/(4*w_mm*1e-3*E_GPa*1e9)) * ((L_mm/h_mm)**3) * array_Fz

plt.errorbar(np.abs(array_Fz), np.abs(array_fleche_mean), array_fleche_std, linestyle='', ecolor='k',color='b',marker='.',label='FEM simulations')
plt.plot(np.abs(array_Fz), np.abs(array_fleche_th),color='tab:orange',label='elastic beam theory')












plt.xlim(0, np.max(np.abs(array_Fz))*1.1)
#plt.ylim(0, np.max(np.abs(array_fleche_mean))*1.1)
plt.xlabel('Fz magnitude (N)')
plt.ylabel('Deflection at the center (m)')
plt.legend()
plt.show()
# %%
# verif que pour 2 modules d'young differents la forme de la plaque déformée est qualitativement la même
# pour ca on normalise les 2 courbes qu'on compare
%matplotlib qt
plt.figure()
plt.plot(dict_results['plate_E3.0GPa_L80mm_w40mm_h4.00mm_Fz-10.0N.vtk']['xmean'], dict_results['plate_E3.0GPa_L80mm_w40mm_h4.00mm_Fz-10.0N.vtk']['uz_mean']/np.max(np.abs(dict_results['plate_E3.0GPa_L80mm_w40mm_h4.00mm_Fz-10.0N.vtk']['uz_mean'])),'.')

plt.plot(dict_results['plate_E9.0GPa_L80mm_w40mm_h4.00mm_Fz-10.0N.vtk']['xmean'], dict_results['plate_E9.0GPa_L80mm_w40mm_h4.00mm_Fz-10.0N.vtk']['uz_mean']/np.max(np.abs(dict_results['plate_E9.0GPa_L80mm_w40mm_h4.00mm_Fz-10.0N.vtk']['uz_mean'])),'+')

plt.plot(dict_results['plate_E3.0GPa_L80mm_w40mm_h2.00mm_Fz-10.0N.vtk']['xmean'], dict_results['plate_E3.0GPa_L80mm_w40mm_h2.00mm_Fz-10.0N.vtk']['uz_mean']/np.max(np.abs(dict_results['plate_E3.0GPa_L80mm_w40mm_h2.00mm_Fz-10.0N.vtk']['uz_mean'])),'+')


plt.show()
# %%
