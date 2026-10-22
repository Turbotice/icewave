#%%
import matplotlib.pyplot as plt
import numpy as np
import time
import cv2

#%%
%matplotlib qt

#%%
# --- Configuration ---
dt = 0.1  # Intervalle de mise à jour (en secondes)
image_source1 = 0  # Caméra 1

# --- Initialisation de la figure avec 2 sous-graphiques ---
fig, ax = plt.subplots()  # 1 ligne, 2 colonnes

# --- Boucle principale ---
t0 = time.time()
while True:
    curr_time = time.time()
    t = curr_time - t0

    # --- Lire les images depuis les 2 caméras ---
    cap1 = cv2.VideoCapture(image_source1)
    ret1, frame1 = cap1.read()
    cap1.release()


    if not ret1:
        print("Erreur : Impossible de lire une des images.")
        break

    # Convertir en RGB
    frame1_rgb = cv2.cvtColor(frame1, cv2.COLOR_BGR2RGB)

    # --- Effacer les anciennes images ---
    for img in ax.images:
        img.remove()

    # --- Afficher les nouvelles images ---
    ax.imshow(frame1_rgb)
    ax.axis('off')
    ax.set_title("Caméra 1")


    # --- Mettre à jour la figure ---
    plt.pause(dt)