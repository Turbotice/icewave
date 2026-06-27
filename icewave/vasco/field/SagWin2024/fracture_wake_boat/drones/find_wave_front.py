import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter

import sys
icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)
from vasco.tools.tridimfits import fit_3dparabola, measure_bidimensional_curvature_around_central_point, convertir_base_tournee, change_order_from_polyfeatures_to_polycoefs



def find_peaks2d(matrice2d, axis=0):
    output_matrix = np.zeros_like(matrice2d) * np.nan
    if axis==0:
        for i in range(matrice2d.shape[1]):
            fp = find_peaks(matrice2d[:,i])[0]
            output_matrix[fp,i] = 1
    elif axis==1:
        for i in range(matrice2d.shape[0]):
            fp = find_peaks(matrice2d[i,:])[0]
            output_matrix[i,fp] = 1
    
    return output_matrix



#matrice2d_exemple = uz[:,:,500]

def find_wave_fronts_on_image(matrice2d, sigma_smooth=5, axis=0, plot=True):

    matrice2d_smoothed = gaussian_filter(matrice2d, sigma_smooth)

    #peaks2d = find_peaks2d(matrice2d_exemple, axis=1)
    peaks2d = find_peaks2d(matrice2d_smoothed, axis=axis)

    if plot:
        plt.figure()
        plt.imshow(matrice2d_smoothed)
        plt.imshow(peaks2d)
        plt.show()

    return peaks2d, matrice2d_smoothed



def find_lines(
    image: np.ndarray,
    min_points_line: int = 2,
    tolerance_px: float = 4.0,
) -> List[List[Tuple[int, int]]]:
    """
    Trouve les lignes formées par les pixels à 1 dans une image binaire.

    Les lignes sont construites par chaînage de plus proches voisins :
    chaque point est relié au point à 1 le plus proche qui n'appartient
    pas encore à une ligne, tant que cette distance reste ≤ tolerance_px.

    Parameters
    ----------
    image : np.ndarray
        Image binaire 2D (valeurs 0 et 1, dtype quelconque).
    min_points_line : int
        Nombre minimal de points pour qu'une séquence soit retournée
        comme une ligne (défaut : 2).
    tolerance_px : float
        Distance euclidienne maximale (en pixels) entre deux points
        consécutifs d'une même ligne (défaut : 4.0).

    Returns
    -------
    List[List[Tuple[int, int]]]
        Liste de lignes ; chaque ligne est une liste ordonnée de points
        (row, col) appartenant à cette ligne.
    """
    if image.ndim != 2:
        raise ValueError("L'image doit être 2D.")

    # Coordonnées (row, col) de tous les pixels à 1
    rows, cols = np.where(image == 1)
    if len(rows) == 0:
        return []

    points = list(zip(rows.tolist(), cols.tolist()))
    available = set(range(len(points)))   # indices des points non encore assignés

    pts_arr = np.array(points, dtype=float)   # shape (N, 2) — pour calculs vectorisés

    lines: List[List[Tuple[int, int]]] = []

    while available:
        # Démarre une nouvelle ligne avec le premier point disponible
        start_idx = next(iter(available))
        available.remove(start_idx)

        current_line = [points[start_idx]]
        current_idx = start_idx

        # Prolonge la ligne par plus proche voisin
        while available:
            current_pt = pts_arr[current_idx]

            # Calcule les distances entre le point courant et tous les disponibles
            avail_idx = np.array(list(available), dtype=int)
            diffs = pts_arr[avail_idx] - current_pt          # (M, 2)
            dists = np.hypot(diffs[:, 0], diffs[:, 1])       # (M,)

            nearest_pos = int(np.argmin(dists))
            nearest_dist = float(dists[nearest_pos])
            nearest_idx = int(avail_idx[nearest_pos])

            if nearest_dist > tolerance_px:
                break   # plus aucun voisin assez proche → fin de cette ligne

            available.remove(nearest_idx)
            current_line.append(points[nearest_idx])
            current_idx = nearest_idx

        if len(current_line) >= min_points_line:
            lines.append(current_line)

    return lines

def group_lines_to_dict(lines=list):
    dict_lines = {}
    for ind_line in range(len(lines)):
        if not(f'line_{ind_line}' in dict_lines):
            dict_lines[f'line_{ind_line}'] = {'x_ind':[], 'y_ind':[]}

        for i in range(len(lines[ind_line])):
            x_ind = lines[ind_line][i][1]
            y_ind = lines[ind_line][i][0]
            dict_lines[f'line_{ind_line}']['x_ind'].append(x_ind)
            dict_lines[f'line_{ind_line}']['y_ind'].append(y_ind)
        
        dict_lines[f'line_{ind_line}']['x_ind'] = np.array(dict_lines[f'line_{ind_line}']['x_ind'])
        dict_lines[f'line_{ind_line}']['y_ind'] = np.array(dict_lines[f'line_{ind_line}']['y_ind'])
        
    
    return dict_lines

#####
#%%
def compute_kappa_n_profiles_along_wavecrest(uz, index_time=1000, facq_x=1.25, plot=True,
                                              params={'sigma_smooth':5, 'axis2dfindpeaks':0,
                                                      'min_points_line':20,'tolerance_px':6.0,
                                                      'fit_higher_order':True}):

    matrice2d = uz[:,:,index_time]

    peaks2d, matrice2d_smoothed = find_wave_fronts_on_image(matrice2d, sigma_smooth=params['sigma_smooth'], axis=params['axis2dfindpeaks'], plot=False)

    if plot:
        plt.figure()
        plt.imshow(matrice2d_smoothed)
        plt.imshow(peaks2d)
        plt.show()


    img = peaks2d

    lines = find_lines(img, min_points_line=params['min_points_line'], tolerance_px=params['tolerance_px'])

    print(f"{len(lines)} ligne(s) trouvée(s) :\n")
    for i, line in enumerate(lines):
        print(f"  Ligne {i} — {len(line)} point(s) : {line}")

    dict_lines = group_lines_to_dict(lines)

    if plot:
        plt.figure()
        plt.imshow(matrice2d, cmap='gray')
    ct=0
    for k in dict_lines:
        xvals2plot = dict_lines[k]['x_ind']
        yvals2plot = dict_lines[k]['y_ind']
        # fit lineaire sur tout le profil
        a,b = np.polyfit(xvals2plot,yvals2plot, 1)

        a4,a3,a2,a1,a0 = np.polyfit(xvals2plot,yvals2plot, 4)
        synth_ydata = a4*xvals2plot**4 + a3*xvals2plot**3 + a2*xvals2plot**2 + a1*xvals2plot + a0

        #synthdata_derivative = 4*a4*xvals2plot**3 + 3*a3*xvals2plot**2 + 2*a2*xvals2plot + a1

        df = np.gradient(synth_ydata)
        dx = np.gradient(xvals2plot)
        
        #arr_alpha = np.arctan2(df, dx)
        #arr_alpha = np.arctan(synthdata_derivative)
        #arr_theta = arr_alpha + np.pi/2

        norm = np.sqrt(1 + (df/dx)**2)
        u_n = -(df/dx) / norm
        v_n = 1/norm
        u_t = 1/norm
        v_t = (df/dx) / norm

        if plot:
            #plt.plot(xvals2plot, a*xvals2plot+b,'k',alpha=0.5)
            plt.plot(xvals2plot, synth_ydata, 'k')
            plt.plot(xvals2plot, yvals2plot, label='line '+str(ct))
            ax = plt.gca()
            if ax.yaxis_inverted():
                varrow = -v_n
            else:
                varrow = v_n
            plt.quiver(
                xvals2plot, yvals2plot, u_n, varrow,
                np.ones(len(xvals2plot)),                    # couleur selon ||∇f||
                cmap='plasma',
                scale=40,               # ajuster la taille des flèches
                width=0.001,
                pivot='mid',
                alpha=0.9
            )

        dict_lines[k]['linear_fit'] = {}
        dict_lines[k]['linear_fit']['slope'] = a
        dict_lines[k]['linear_fit']['intercept'] = b

        dict_lines[k]['polynomial_fit'] = {}
        dict_lines[k]['polynomial_fit']['a4'] = a4
        dict_lines[k]['polynomial_fit']['a3'] = a3
        dict_lines[k]['polynomial_fit']['a2'] = a2
        dict_lines[k]['polynomial_fit']['a1'] = a1
        dict_lines[k]['polynomial_fit']['a0'] = a0

        dict_lines[k]['polynomial_fit']['xvalsfit'] = xvals2plot
        dict_lines[k]['polynomial_fit']['yfit'] = synth_ydata

        dict_lines[k]['polynomial_fit']['u_n'] = u_n
        dict_lines[k]['polynomial_fit']['v_n'] = v_n
        dict_lines[k]['polynomial_fit']['u_t'] = u_t
        dict_lines[k]['polynomial_fit']['v_t'] = v_t
        


        ct+=1
    if plot:
        plt.legend()
        plt.show()


    # tracé de la courbure "transverse" le long d'une crête de vague

    # tout d'abord, on peut calculer le champ de courbures :

    #kappa_x = np.gradient(np.gradient(matrice2d_smoothed,axis=1),axis=1) * (facq_x**2)
    #kappa_y = np.gradient(np.gradient(matrice2d_smoothed,axis=0),axis=0) * (facq_x**2)

    kappa_x = np.diff(np.diff(matrice2d_smoothed,axis=1,prepend=0),axis=1,prepend=0) * (facq_x**2)
    kappa_y = np.diff(np.diff(matrice2d_smoothed,axis=0,prepend=0),axis=0,prepend=0) * (facq_x**2)

    
    if plot:
        plt.figure()
    for k in dict_lines:

        x_ind = dict_lines[k]['x_ind']
        y_ind = dict_lines[k]['y_ind']
        kappa_n_prof = np.zeros(len(x_ind))
        kappa_t_prof = np.zeros(len(x_ind))
        kappa_x_prof = np.zeros(len(x_ind))
        kappa_y_prof = np.zeros(len(x_ind))

        if not(params['fit_higher_order']):

            slope = dict_lines[k]['linear_fit']['slope']
            theta = np.arctan(slope) + np.pi/2
            alpha = np.arctan(slope)
            dict_lines[k]['linear_fit']['theta_rad'] = theta
            for i in range(len(x_ind)):
                kappa_x_prof[i] = kappa_x[y_ind[i],x_ind[i]]
                kappa_y_prof[i] = kappa_y[y_ind[i],x_ind[i]]
                kappa_n_prof[i] = np.cos(theta) * kappa_x[y_ind[i],x_ind[i]] + np.sin(theta) * kappa_y[y_ind[i],x_ind[i]]
                kappa_t_prof[i] = np.cos(alpha) * kappa_x[y_ind[i],x_ind[i]] + np.sin(alpha) * kappa_y[y_ind[i],x_ind[i]]

        elif params['fit_higher_order']:
            u_n = dict_lines[k]['polynomial_fit']['u_n']
            v_n = dict_lines[k]['polynomial_fit']['v_n']
            u_t = dict_lines[k]['polynomial_fit']['u_t']
            v_t = dict_lines[k]['polynomial_fit']['v_t']
            

            for i in range(len(x_ind)):
                kappa_x_prof[i] = kappa_x[y_ind[i],x_ind[i]]
                kappa_y_prof[i] = kappa_y[y_ind[i],x_ind[i]]
                kappa_n_prof[i] = u_n[i] * kappa_x[y_ind[i],x_ind[i]] + v_n[i] * kappa_y[y_ind[i],x_ind[i]]
                kappa_t_prof[i] = u_t[i] * kappa_x[y_ind[i],x_ind[i]] + v_t[i] * kappa_y[y_ind[i],x_ind[i]]


        dict_lines[k]['kappa_n_profile'] = kappa_n_prof
        dict_lines[k]['kappa_t_profile'] = kappa_t_prof
        dict_lines[k]['kappa_x_profile'] = kappa_x_prof
        dict_lines[k]['kappa_y_profile'] = kappa_y_prof
        
        
        if plot:  
            plt.plot(-kappa_n_prof)
            plt.ylim(0,None)
            

    return dict_lines, (peaks2d, matrice2d_smoothed, matrice2d)

#%%
def bidim_curvature_oneline(uz=np.ndarray,
                            dic_all_lines=dict,
                            ind_t=1088, 
                            line_index=0,
                            window_size=(11,11)):
    """
    info : pour ind_t il faut mettre ind_tfrac_approx correspondant à 
    l'instant approx d'une fracture
    """
    
    xind_arr = dic_all_lines[f'frame_{ind_t}']['dict_lines'][f'line_{line_index}']['x_ind']
    yind_arr = dic_all_lines[f'frame_{ind_t}']['dict_lines'][f'line_{line_index}']['y_ind']
    u_n_arr = dic_all_lines[f'frame_{ind_t}']['dict_lines'][f'line_{line_index}']['polynomial_fit']['u_n']
    u_t_arr = dic_all_lines[f'frame_{ind_t}']['dict_lines'][f'line_{line_index}']['polynomial_fit']['u_t']
    v_n_arr = dic_all_lines[f'frame_{ind_t}']['dict_lines'][f'line_{line_index}']['polynomial_fit']['v_n']
    v_t_arr = dic_all_lines[f'frame_{ind_t}']['dict_lines'][f'line_{line_index}']['polynomial_fit']['v_t']

    

    alpha_arr = np.zeros(len(u_t_arr))
    for i in range(len(u_t_arr)):
        alpha_arr[i] = np.arctan2(v_t_arr[i], u_t_arr[i])
    print(alpha_arr)
    
    arr_coefficients_invariantbasis = np.zeros((len(xind_arr), 6)) # 6 coef pour un fit de parabole tridimensionnelle
    arr_coefficients_movingbasis = np.zeros((len(xind_arr), 6))

    for index_point in range(len(xind_arr)):
        if ((xind_arr[index_point]-int((window_size[0]-1)/2) < 0) | (yind_arr[index_point]-int((window_size[1]-1)/2) < 0)):
            arr_coefficients_invariantbasis[i,:] = np.ones(6)*np.nan
            arr_coefficients_movingbasis[i,:] = np.ones(6)*np.nan
        else:
            matrix2d = uz[:,:,ind_t]
            dic_res = measure_bidimensional_curvature_around_central_point(matrix2d, xind_arr[index_point], yind_arr[index_point], window_size=window_size)
            coef_polyfeatures = dic_res['coefficients']
            intercept = dic_res['intercept']

            coefficients_invariantbasis = change_order_from_polyfeatures_to_polycoefs(coef_polyfeatures=coef_polyfeatures, intercept=intercept)
            arr_coefficients_invariantbasis[index_point,:] = coefficients_invariantbasis

            coefficients_movingbasis = convertir_base_tournee(coefficients_invariantbasis, alpha_arr[index_point])
            arr_coefficients_movingbasis[index_point,:] = coefficients_movingbasis

    return arr_coefficients_invariantbasis, arr_coefficients_movingbasis, alpha_arr