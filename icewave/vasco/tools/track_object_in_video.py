#%%
import cv2
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
import pickle
#%%
def load_image_sequence(folder_path, extension="tiff"):
    pattern = os.path.join(folder_path, f"*.{extension}")
    files = sorted(glob.glob(pattern))

    if len(files) == 0:
        raise ValueError("Aucune image trouvée dans le dossier")

    return files

# =========================
# 1. Rotation d'image
# =========================

def rotate_image(image, angle):
    h, w = image.shape
    center = (w // 2, h // 2)

    M = cv2.getRotationMatrix2D(center, angle, 1.0)
    rotated = cv2.warpAffine(image, M, (w, h))

    return rotated


# =========================
# 2. Corrélation 2D
# =========================
"""def compute_correlation(image, template):
    result = cv2.matchTemplate(image, template, cv2.TM_CCORR_NORMED)
    return result


def extract_max_location_subpix(corr_map, debug=True):
    y, x = np.unravel_index(np.argmax(corr_map), corr_map.shape)

    def refine(vec, i):
        if i < 1 or i >= len(vec) - 1:
            return float(i)
        a, b, c = np.log(vec[i-1]), np.log(vec[i]), np.log(vec[i+1])
        denom = 2 * (a + c - 2*b)
        return i + (a - c) / denom if denom else float(i)

    x_sub = refine(corr_map[y, :], x)
    y_sub = refine(corr_map[:, x], y)

    if debug:
        plt.figure(figsize=(6, 5))
        plt.imshow(corr_map, cmap='hot', origin='upper')
        plt.colorbar(label='corrélation')
        plt.scatter(x, y, c='cyan', s=80, marker='+', linewidths=2, label=f'max pixel ({x}, {y})')
        plt.scatter(x_sub, y_sub, c='lime', s=80, marker='x', linewidths=2, label=f'subpixel ({x_sub:.2f}, {y_sub:.2f})')
        plt.legend(fontsize=8)
        plt.title('Carte de corrélation')
        plt.tight_layout()
        plt.show()

    return (y_sub, x_sub), corr_map[y, x]"""

def compute_correlation(image, template):
    result = cv2.matchTemplate(image, template, cv2.TM_CCORR_NORMED)
    return result


def extract_max_location_subpix(corr_map, debug=False):
    y, x = np.unravel_index(np.argmax(corr_map), corr_map.shape)

    def refine_gaussian(vec, i):
        if i < 1 or i >= len(vec) - 1:
            return float(i), vec[i], None
        a, b, c = np.log(vec[i-1]), np.log(vec[i]), np.log(vec[i+1])
        denom = 2 * (a + c - 2*b)
        if not denom:
            return float(i), vec[i], None
        delta = (a - c) / denom
        i_sub = i + delta
        log_A = b - (a - c)**2 / (4 * denom)  # = b - delta*(a-c)/4
        sigma2 = -1 / denom
        return i_sub, log_A, sigma2

    x_sub, logAx, sigma2x = refine_gaussian(corr_map[y, :], x)
    y_sub, logAy, sigma2y = refine_gaussian(corr_map[:, x], y)

    b = np.log(corr_map[y, x])

    # A est compté deux fois via b dans logAx et logAy
    log_A = logAx + logAy - b
    interp_val = np.exp(log_A)

    if debug:
        plt.figure(figsize=(6, 5))
        plt.imshow(corr_map, cmap='hot', origin='upper')
        plt.colorbar(label='corrélation')
        plt.scatter(x, y, c='cyan', s=80, marker='+', linewidths=2, label=f'max pixel ({x}, {y})')
        plt.scatter(x_sub, y_sub, c='lime', s=80, marker='x', linewidths=2, label=f'subpixel ({x_sub:.2f}, {y_sub:.2f})')
        plt.legend(fontsize=8)
        plt.title('Carte de corrélation')
        plt.tight_layout()
        plt.show()

    return (y_sub, x_sub), interp_val

def extract_max_location(corr_map):
    (y_sub, x_sub), max_val = extract_max_location_subpix(corr_map)
    return max_val, (int(x_sub), int(y_sub))

# =========================
# 4. Recherche sur rotation
# =========================
def search_best_transform(frame, template, angles):
    best_score = -np.inf
    best_params = (0, 0, 0)

    for angle in angles:
        #rotated_template = rotate_image(template, angle)
        rotated_frame = rotate_image(frame, -angle)

        #corr_map = compute_correlation(frame, rotated_template)
        corr_map = compute_correlation(rotated_frame, template)
        (y, x), score = extract_max_location_subpix(corr_map)

        if score > best_score:
            best_score = score
            best_params = (x, y, angle)

    return best_params, best_score


# =========================
# 5. Mise à jour template
# =========================
"""def extract_new_template(frame, position, size):
    x, y = position
    w, h = size

    return frame[y:y+h, x:x+w]
"""
def extract_new_template(frame, position, size):
    x, y = position  # floats
    w, h = size

    # Matrice de translation subpixel
    M = np.float32([
        [1, 0, -x],
        [0, 1, -y]
    ])

    # Extraction avec interpolation bilinéaire
    patch = cv2.warpAffine(
        frame, M, (w, h),
        flags=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_REFLECT
    )

    return patch

def resize_for_display(frame, max_width=1700, max_height=800):
    h, w = frame.shape[:2]

    scale = min(max_width / w, max_height / h)

    if scale < 1:
        new_size = (int(w * scale), int(h * scale))
        return cv2.resize(frame, new_size)
    return frame

def frame_generator(video_path=None, image_files=None):

    if video_path is not None:
        cap = cv2.VideoCapture(video_path)

        while True:
            ret, frame = cap.read()
            if not ret:
                break
            yield frame

        cap.release()

    elif image_files is not None:
        for file in image_files:
            frame = cv2.imread(file)
            if frame is None:
                continue
            yield frame

    else:
        raise ValueError("Fournir soit video_path soit image_files")

# =========================
# 6. Tracker principal
# =========================
def track_object(video_path=None, image_files=None, initial_bbox=(0,0,0,0), angle_range=(-10,10,21)):

    frames = frame_generator(video_path, image_files)

    # première frame
    frame = next(frames, None)
    if frame is None:
        print("Erreur lecture")
        return

    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

    x, y, w, h = initial_bbox
    template = gray[y:y+h, x:x+w]

    angles = np.linspace(angle_range[0], angle_range[1], angle_range[2])
    results = []
    count = 0
    for frame in frames:
        print('iteration',count)
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        (dx, dy, theta), score = search_best_transform(gray, template, angles)

        results.append({
            "x": dx,
            "y": dy,
            "theta": theta,
            "score": score
        })
        top_left = (np.round(dx).astype(int), np.round(dy).astype(int))
        bottom_right = (np.round(dx + w).astype(int), np.round(dy + h).astype(int))

        display = frame.copy()
        cv2.rectangle(display, top_left, bottom_right, (0, 255, 0), 2)

        cv2.putText(display, f"theta={theta:.2f}", (np.round(dx).astype(int), np.round(dy).astype(int) - 10),
                    cv2.FONT_HERSHEY_SIMPLEX, 0.6, (0,255,0), 2)

        display_resized = resize_for_display(display)
        cv2.imshow("Tracking", display_resized)

        template = extract_new_template(gray, (dx, dy), (w, h))

        if cv2.waitKey(30) & 0xFF == 27:
            break
        count+=1

    cv2.destroyAllWindows()
    return results

def plot_variable(results, var_name="x", fps=1):
    """
    Trace une variable (x, y, theta, score...) en fonction du temps
    
    results : liste de dictionnaires
    var_name : nom de la variable à tracer
    fps : images par seconde
    """

    if len(results) == 0:
        print("Aucun résultat à afficher")
        return

    if var_name not in results[0]:
        print(f"Variable '{var_name}' non trouvée")
        return

    values = [r[var_name] for r in results]
    t = [i / fps for i in range(len(values))]

    plt.figure()
    plt.plot(t, values)
    if fps==1:
        plt.xlabel('Time [frames]')
    else:
        plt.xlabel("Temps (s)")
    plt.ylabel(var_name)
    plt.title(f"Evolution de {var_name} dans le temps")
    plt.grid()

    plt.show()

    return t, values

#%%
if __name__ == "__main__":
    saveresults = True
    plotresults = True
    gendir = "C:/Users/Vasco Zanchi/Desktop/piv_brazilian_test/P1047763"
    base = f"{gendir}/image_sequence"  # dossier contenant les images
    image_files = load_image_sequence(base, "jpg")

    initial_bbox = (1030, 110, 1840, 1840)
    angle_range = (-1, 1, 21)

    results = track_object(
        image_files=image_files,
        initial_bbox=initial_bbox,
        angle_range=angle_range
    )
    if plotresults:
        plot_variable(results, "x")
        plot_variable(results, "y")
        plot_variable(results, "theta")
        plot_variable(results, "score")
    if saveresults:
        output_file = f"{gendir}/results.pkl"
        with open(output_file, 'wb') as file:
            pickle.dump(results, file)
