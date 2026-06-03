#%%
import numpy as np
import matplotlib.pyplot as plt

from skimage.measure import profile_line
import sys

icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)
from vasco.tools import clickonfigures

#%%
import numpy as np
from skimage.measure import profile_line


def all_profiles(image, angle_deg, step=1, linewidth=1, order=1):
    """
    Calcule tous les profils parallèles traversant une image.

    Parameters
    ----------
    image : ndarray (H, W)
        Image 2D.
    angle_deg : float
        Angle des profils en degrés (0° = horizontal vers la droite).
    step : float, optional
        Espacement entre profils en pixels. Défaut = 1.
    linewidth : int, optional
        Paramètre transmis à profile_line.
    order : int, optional
        Ordre d'interpolation transmis à profile_line.

    Returns
    -------
    profiles : list of dict
        Chaque élément contient :
            - 'profile' : valeurs du profil
            - 'src' : point de départ (row, col)
            - 'dst' : point d'arrivée (row, col)
            - 'offset' : décalage perpendiculaire
    """

    h, w = image.shape[:2]

    theta = np.deg2rad(angle_deg)

    # Direction du profil
    d = np.array([np.sin(theta), np.cos(theta)])

    # Direction perpendiculaire
    n = np.array([-d[1], d[0]])

    # Centre de l'image
    center = np.array([(h - 1) / 2, (w - 1) / 2])

    corners = np.array([
        [0, 0],
        [0, w - 1],
        [h - 1, 0],
        [h - 1, w - 1]
    ])

    # Étendue maximale du balayage
    offsets = np.dot(corners - center, n)
    offset_min = offsets.min()
    offset_max = offsets.max()

    profiles = []

    for offset in np.arange(offset_min, offset_max + step, step):

        p0 = center + offset * n

        intersections = []

        # Intersections avec les bords horizontaux
        if abs(d[0]) > 1e-12:
            for r in [0, h - 1]:
                t = (r - p0[0]) / d[0]
                c = p0[1] + t * d[1]
                if 0 <= c <= w - 1:
                    intersections.append([r, c])

        # Intersections avec les bords verticaux
        if abs(d[1]) > 1e-12:
            for c in [0, w - 1]:
                t = (c - p0[1]) / d[1]
                r = p0[0] + t * d[0]
                if 0 <= r <= h - 1:
                    intersections.append([r, c])

        # Suppression des doublons
        unique = []
        for p in intersections:
            if not any(np.allclose(p, q, atol=1e-8) for q in unique):
                unique.append(p)

        if len(unique) < 2:
            continue

        pts = np.array(unique)

        # Paramètre le long de la ligne
        tvals = np.dot(pts - p0, d)

        src = pts[np.argmin(tvals)]
        dst = pts[np.argmax(tvals)]

        # Vérifie que la ligne couvre au moins 2 pixels
        length = np.linalg.norm(dst - src)
        if length < 1:
            continue

        prof = profile_line(
            image,
            src,
            dst,
            linewidth=linewidth,
            order=order,
            mode='constant'
        )

        profiles.append({
            "profile": prof,
            "src": tuple(src),
            "dst": tuple(dst),
            "offset": offset
        })

    return profiles
# %%

def plot_profilelines_onimage(image, profiles):
    plt.figure(figsize=(8, 8))
    plt.imshow(image, cmap="gray")

    for p in profiles:
        src = p["src"]  # (row, col)
        dst = p["dst"]

        plt.plot(
            [src[1], dst[1]],  # x = colonnes
            [src[0], dst[0]],  # y = lignes
            'r-',
            linewidth=0.5,
            alpha=0.3
        )

    plt.axis("image")
    plt.xlim(0,image.shape[1])
    plt.ylim(0,image.shape[0])
    plt.show()
    
"""
profiles = all_profiles(image=image, angle_deg=-40)
plot_profilelines_onimage(image,profiles)"""