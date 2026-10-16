#%%
import os

import cv2
import numpy as np

from PIL import Image
import matplotlib.pyplot as plt

#%%
def open_image(image_path):
    # Ouvrir l'image avec Pillow (pour gérer le TIFF)
    img_pil = Image.open(image_path)

    # Convertir en tableau NumPy et en niveaux de gris
    img_cv = np.array(img_pil)
    if len(img_cv.shape) == 3:  # Si l'image est en couleur (BGR ou RGB)
        img_gray = cv2.cvtColor(img_cv, cv2.COLOR_RGB2GRAY)
    else:
        img_gray = img_cv  # Déjà en niveaux de gris
    
    return img_gray

def appliquer_seuil(img_gray, seuil=100, max_val=255, type_seuil=cv2.THRESH_BINARY):
    """
    Applique un seuil à une image

    Args:
        Image (matrice numpy.ndarray)
        seuil (int): Valeur du seuil (0-255). Par défaut : 127.
        max_val (int): Valeur maximale pour les pixels au-dessus du seuil. Par défaut : 255.
        type_seuil: Type de seuillage OpenCV. Options :
            - cv2.THRESH_BINARY (par défaut)
            - cv2.THRESH_BINARY_INV
            - cv2.THRESH_TRUNC
            - cv2.THRESH_TOZERO
            - cv2.THRESH_TOZERO_INV

    Returns:
        numpy.ndarray: Image seuillée en niveaux de gris (1 canal).
    """

    # Appliquer le seuil
    _, img_seuillee = cv2.threshold(
        img_gray,
        seuil,
        max_val,
        type_seuil
    )

    return img_seuillee

def crop_image(img,xinf,xsup,yinf,ysup):
    return img[yinf:ysup,xinf:xsup] 

def mask_image_between4pts(img, pt1, pt2, pt3, pt4, fill_value=0):
    """
    Met à `fill_value` tous les pixels en dehors du quadrilatère.
    Retourne une image de même type que l'entrée (uint8 par défaut).

    Args:
        img: Image d'entrée (2D ou 3D).
        pt1, pt2, pt3, pt4: Points (x, y) du quadrilatère.
        fill_value: Valeur à appliquer à l'extérieur (0 par défaut).
    """
    h, w = img.shape[:2]
    mask = np.zeros((h, w), dtype=np.uint8)

    # Remplir le quadrilatère en blanc (255)
    pts = np.array([pt1, pt2, pt3, pt4], dtype=np.int32)
    cv2.fillPoly(mask, [pts], color=255)

    # Créer une copie de l'image et appliquer le masque
    masked_img = img.copy()
    if len(img.shape) == 2:  # Niveaux de gris
        masked_img[mask == 0] = fill_value
    else:  # Couleur
        masked_img[mask == 0] = fill_value  # fill_value doit être un tuple pour BGR (ex: (0,0,0))

    return masked_img



def rotate_bound(image, angle):
    """
    Tourne une image sans couper ses bords.
    La taille de l'image résultante est adaptée à la rotation.
    """

    h, w = image.shape[:2]

    center = (w / 2, h / 2)

    # Matrice de rotation
    M = cv2.getRotationMatrix2D(
        center,
        angle,
        1.0
    )

    # Calcul de la nouvelle taille
    cos = abs(M[0, 0])
    sin = abs(M[0, 1])

    new_w = int((h * sin) + (w * cos))
    new_h = int((h * cos) + (w * sin))

    # Ajuster la translation pour centrer l'image
    M[0, 2] += (new_w / 2) - center[0]
    M[1, 2] += (new_h / 2) - center[1]

    # Rotation
    rotated = cv2.warpAffine(
        image,
        M,
        (new_w, new_h),
        flags=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=0
    )

    return rotated

def match_with_rotation(
    roi,
    template,
    angle_range=10,
    angle_step=1,
    coarse_step=3,
    show_correlation=False,
    show_match=False,
    mask_background=False,
    black_weight=0.1
):
    """
    Recherche un template dans une ROI en deux étapes :
    recherche grossière, puis affinage angulaire.
    """

    if roi.ndim != 2 or template.ndim != 2:
        raise ValueError("roi et template doivent être en grayscale 2D.")

    rh, rw = roi.shape

    if template.shape[0] > rh or template.shape[1] > rw:
        raise ValueError("Le template doit être plus petit que la ROI.")

    best_score = -np.inf
    best_angle = 0
    best_pos = None
    best_size = None
    best_corr_map = None
    best_template = None

    # Cache les résultats pour ne pas recalculer un angle déjà testé.
    tested = {}

    def test_angle(angle):
        nonlocal best_score, best_angle, best_pos
        nonlocal best_size, best_corr_map, best_template

        # Un angle déjà évalué n'est pas recalculé.
        key = round(float(angle), 6)
        if key in tested:
            return

        rotated = rotate_bound(template, angle)
        th, tw = rotated.shape[:2]

        if th > rh or tw > rw:
            tested[key] = None
            return

        if mask_background:
            mask = np.where(
                rotated > 127,
                1.0,
                black_weight
            ).astype(np.float32)

            corr = cv2.matchTemplate(
                roi, rotated, cv2.TM_CCORR_NORMED, mask=mask
            )
        else:
            corr = cv2.matchTemplate(
                roi, rotated, cv2.TM_CCOEFF_NORMED
            )

        valid = np.isfinite(corr)
        if not np.any(valid):
            tested[key] = None
            return

        safe_corr = np.where(valid, corr, -np.inf)
        _, score, _, pos = cv2.minMaxLoc(safe_corr)

        tested[key] = (score, pos, (tw, th), corr, rotated)

        if score > best_score:
            best_score = score
            best_angle = float(angle)
            best_pos = pos
            best_size = (tw, th)
            best_corr_map = corr
            best_template = rotated.copy()

    # Étape 1 : recherche grossière
    coarse_step = max(float(coarse_step), float(angle_step))
    coarse_angles = np.arange(
        -angle_range,
        angle_range + coarse_step * 0.5,
        coarse_step
    )

    for angle in coarse_angles:
        test_angle(angle)

    if best_pos is None:
        return -1, 0, None, None

    # Étape 2 : affiner autour du meilleur angle grossier
    coarse_best_angle = best_angle
    fine_angles = np.arange(
        max(-angle_range, coarse_best_angle - coarse_step),
        min(angle_range, coarse_best_angle + coarse_step) + angle_step * 0.5,
        angle_step
    )

    for angle in fine_angles:
        test_angle(angle)

    # Affichage de la carte de corrélation
    if show_correlation and best_corr_map is not None:
        plt.figure(figsize=(10, 6))
        plt.imshow(best_corr_map, cmap="viridis")
        plt.colorbar(label="Corrélation")
        plt.title(f"Score={best_score:.3f}, angle={best_angle:.1f}°")
        plt.xlabel("X")
        plt.ylabel("Y")
        plt.show()

    # Affichage du meilleur match
    if show_match and best_template is not None:
        x, y = best_pos
        w, h = best_size

        plt.figure(figsize=(10, 6))
        plt.imshow(roi, cmap="gray")
        plt.imshow(
            best_template,
            cmap="jet",
            alpha=0.6,
            extent=[x, x + w, y + h, y]
        )

        rectangle = plt.Rectangle(
            (x, y), w, h,
            fill=False, edgecolor="red", linewidth=2
        )
        plt.gca().add_patch(rectangle)
        plt.xlim(0, rw)
        plt.ylim(rh, 0)
        plt.title(
            f"Score={best_score:.3f}, angle={best_angle:.1f}°, "
            f"position=({x}, {y})"
        )
        plt.show()

    return best_score, best_angle, best_pos, best_size


def create_digit_templates(
    size=(520, 190),
    thickness=50,
    color=255,
    padding=0
):
    """
    Crée des templates de digits 7 segments avec une ROI plus grande
    contenant des zéros autour du digit.

    Args:
        size: (hauteur, largeur) du digit.
        thickness: épaisseur des segments.
        color: couleur des segments.
        padding: nombre de pixels ajoutés autour du digit.

    Returns:
        templates: dictionnaire contenant les templates 0-9.
    """

    templates = {}

    # Taille du digit
    h, w = size

    # Taille totale de la ROI
    roi_h = h + 2 * padding
    roi_w = w + 2 * padding

    # Définition des segments pour chaque digit
    digit_segments = {
        0: [1, 1, 1, 1, 1, 1, 0],
        1: [0, 1, 1, 0, 0, 0, 0],
        2: [1, 1, 0, 1, 1, 0, 1],
        3: [1, 1, 1, 1, 0, 0, 1],
        4: [0, 1, 1, 0, 0, 1, 1],
        5: [1, 0, 1, 1, 0, 1, 1],
        6: [1, 0, 1, 1, 1, 1, 1],
        7: [1, 1, 1, 0, 0, 0, 0],
        8: [1, 1, 1, 1, 1, 1, 1],
        9: [1, 1, 1, 1, 0, 1, 1],
    }

    # Coordonnées des segments
    segment_coords = {
        'a': [
            (w // 10, h // 10),
            (9 * w // 10, h // 10)
        ],

        'b': [
            (9 * w // 10, h // 10),
            (9 * w // 10, h // 2)
        ],

        'c': [
            (9 * w // 10, h // 2),
            (9 * w // 10, 9 * h // 10)
        ],

        'd': [
            (w // 10, 9 * h // 10),
            (9 * w // 10, 9 * h // 10)
        ],

        'e': [
            (w // 10, h // 2),
            (w // 10, 9 * h // 10)
        ],

        'f': [
            (w // 10, h // 10),
            (w // 10, h // 2)
        ],

        'g': [
            (w // 10, h // 2),
            (9 * w // 10, h // 2)
        ],
    }

    # Créer chaque digit
    for digit, segments in digit_segments.items():

        # Décalage horizontal spécifique au digit 1
        offset_x = (9 * w // 10) - (w // 2) if digit == 1 else 0

        # Grande ROI noire
        template = np.zeros(
            (roi_h, roi_w),
            dtype=np.uint8
        )

        # Dessiner les segments avec un décalage
        # correspondant au padding
        for seg_name, is_on in zip(
            ['a', 'b', 'c', 'd', 'e', 'f', 'g'],
            segments
        ):

            if is_on:

                start, end = segment_coords[seg_name]

                # Décalage dans la grande ROI
                start = (
                    int(start[0] + padding - offset_x),
                    int(start[1] + padding)
                )

                end = (
                    int(end[0] + padding - offset_x),
                    int(end[1] + padding)
                )

                cv2.line(
                    template,
                    start,
                    end,
                    color,
                    thickness
                )

        templates[digit] = template

    return templates


def detect_multiple_digits(
    image,
    templates,
    threshold=0.60,
    angle_range=10,
    angle_step=2,
    nms_threshold=0.3,
    mask_background=False,
    black_weight=0.1
):
    """
    Détecte plusieurs occurrences de chiffres dans une image.

    Retourne une liste de détections :
    digit, score, angle, position et taille.
    """

    detections = []

    for digit, template in templates.items():
        for angle in np.arange(
            -angle_range,
            angle_range + angle_step * 0.5,
            angle_step
        ):
            rotated = rotate_bound(template, angle)
            th, tw = rotated.shape[:2]

            if th > image.shape[0] or tw > image.shape[1]:
                continue

            if mask_background:
                mask = np.where(
                    rotated > 127, 1.0, black_weight
                ).astype(np.float32)

                corr = cv2.matchTemplate(
                    image,
                    rotated,
                    cv2.TM_CCORR_NORMED,
                    mask=mask
                )
            else:
                corr = cv2.matchTemplate(
                    image,
                    rotated,
                    cv2.TM_CCOEFF_NORMED
                )

            # Remplacer les valeurs invalides
            corr = np.where(np.isfinite(corr), corr, -np.inf)

            # Maxima locaux : pixels qui sont les plus élevés
            # dans leur voisinage.
            local_max = cv2.dilate(
                corr.astype(np.float32),
                np.ones((th // 2 or 1, tw // 2 or 1), np.uint8)
            )

            ys, xs = np.where(
                (corr >= threshold) & (corr == local_max)
            )

            for x, y in zip(xs, ys):
                detections.append({
                    "digit": digit,
                    "score": float(corr[y, x]),
                    "angle": float(angle),
                    "position": (int(x), int(y)),
                    "size": (tw, th)
                })

    # Tri : meilleurs scores en premier
    detections.sort(key=lambda d: d["score"], reverse=True)

    # Suppression des détections qui se chevauchent trop (NMS)
    kept = []

    for det in detections:
        x, y = det["position"]
        w, h = det["size"]
        area = w * h
        duplicate = False

        for prev in kept:
            px, py = prev["position"]
            pw, ph = prev["size"]

            ix1 = max(x, px)
            iy1 = max(y, py)
            ix2 = min(x + w, px + pw)
            iy2 = min(y + h, py + ph)

            intersection = max(0, ix2 - ix1) * max(0, iy2 - iy1)
            union = area + pw * ph - intersection

            iou = intersection / union if union else 0

            if iou > nms_threshold:
                duplicate = True
                break

        if not duplicate:
            kept.append(det)

    return kept

def creer_sequence(detections, top_n=3):
    # Garder les top_n détections ayant les meilleurs scores
    meilleures = sorted(
        detections,
        key=lambda d: d["score"],
        reverse=True
    )[:top_n]

    # Trier ces détections de gauche à droite selon x
    meilleures = sorted(
        meilleures,
        key=lambda d: d["position"][0]
    )

    # Construire la séquence de chiffres
    sequence = "".join(str(d["digit"]) for d in meilleures)

    return sequence

def preprocess_image(image,
                     seuil=100,
                     corners_mask = [(816,1415), (814,945), (1789,1011), (1694,1539)]):
    image_seuillee = appliquer_seuil(image, seuil=100)
    image_seuillee_inversee = cv2.bitwise_not(image_seuillee)

    %matplotlib inline


    if type(corners_mask)==list:
        image_preprocessed = mask_image_between4pts(image_seuillee_inversee,corners_mask[0],corners_mask[1],corners_mask[2],corners_mask[3])
    else:
        image_preprocessed = image_seuillee_inversee

    return image_preprocessed


#%%

if __name__=='__main__':
    # test de la detection de digits

    templates = create_digit_templates()

    img_path = r"C:\Users\Vasco Zanchi\Desktop\CASCADE26\capteur_deplacement\test1_facq19Hz_expo5000\Basler_acA3088-57um__23965290__20261014_163815303_0085.tiff"
    
    image = open_image(image_path=img_path)
    corners_mask = [(816,1415), (814,945), (1789,1011), (1694,1539)]
    image_preprocessed = preprocess_image(image=image,seuil=100, corners_mask=corners_mask)
    
    detections = detect_multiple_digits(
            cv2.bitwise_not(image_preprocessed),
            templates,
            threshold=0.6,
            angle_range=10,
            angle_step=2
        )

    for det in detections:
        print(
            det["digit"],
            det["score"],
            det["position"],
            det["angle"]
        )


    # test :
    sequence = creer_sequence(detections, top_n=3)
    print("Séquence détectée :", sequence)

    plt.imshow(image_preprocessed)

    for i in range(len(detections)):
        if i<3:
            plt.plot(detections[i]['position'][0], detections[i]['position'][1], 'ro')
        else:
            break
    plt.title('Detected digits sequence : ' + sequence)
# %%
