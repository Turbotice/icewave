import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple

from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter

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