#%% imports modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

import pickle
import os
import sys

icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

import tools.rw_data as rw

from vasco.tools.clickonfigures import profile_line_on_image_2clicks
from vasco.tools.clickonfigures import get_n_points_onimage
from vasco.tools.clickonfigures import get_n_points

from functions_fracture_analysis import click_on_fracture_path_plot_time_evol, click2extract_amplitude, plot_elevation_refnotbroken_and_broken
from profiles import *

#%% definition des chemins des données
disk = 'L:'# disk is Elements on adour
date = '0211'

data_path = f'{disk}/Share_hublot/Data'
daily_drone_data_path = f'{data_path}/{date}/Drones'

velocity_field_path = f'{daily_drone_data_path}/exact_solution_real_field_stereo_0211_2024_rectangular_grid.h5'

# traitement stéphane : postitions des fractures (obtenue avec divergence du champ de vitesses)
fractures_positions_path = f'{daily_drone_data_path}/Results/fracture_positions.pkl'

#%% Cellule à exécuter une ceule fois (chargement champs vitesses from stereo piv)
#load dict from h5 file

dict_stereo_pivdata = rw.load_dict_from_h5(velocity_field_path)

#%% chargement des positions des fractures détectées
with open(fractures_positions_path, 'rb') as file:
    fractures_positions_data = pickle.load(file)

# %% definitions des variables utiles (champs d'élévation etc.)
dict_stereo_pivdata.keys()

vx = dict_stereo_pivdata['u'][0,:,:,:]
vy = dict_stereo_pivdata['u'][1,:,:,:]
vz = dict_stereo_pivdata['u'][2,:,:,:]

dt = dict_stereo_pivdata['t'][1] - dict_stereo_pivdata['t'][0]

ux = np.cumsum(vx, axis=2)*dt
uy = np.cumsum(vy, axis=2)*dt
uz = np.cumsum(vz, axis=2)*dt

facq_x = dict_stereo_pivdata['SCALE']['facq_x']

#%%

# à continuer pour inclure les positions des fractures etc.


#%%

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


output_video_path = f'{daily_drone_data_path}/Results/traitement_vasco/video_uz_test.mp4'



# ── Vos données : tableau numpy (Ny, Nx, Nt) ──────────────────────────────────

y = np.arange(uz.shape[0])
x = np.arange(uz.shape[1])
t = np.arange(uz.shape[2])
field = uz

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 5.6))
fig.subplots_adjust(right=0.88)

vmin, vmax = field.min(), field.max()

im = ax.imshow(
    field[:, :, 0],
    origin="lower",
    cmap="RdBu_r",
    vmin=vmin, vmax=vmax,
    aspect="auto",
    extent=[x.min(), x.max(), y.min(), y.max()],
)
fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("Amplitude")
ax.set_xlabel("x")
ax.set_ylabel("y")
title = ax.set_title("t = 0.00")

# ── Update ────────────────────────────────────────────────────────────────────
def update(frame):
    im.set_data(field[:, :, frame])
    title.set_text(f"t = {t[frame]:.2f}")
    return im, title

# ── Animation ─────────────────────────────────────────────────────────────────
ani = animation.FuncAnimation(fig, update, frames=len(t), interval=1000/20, blit=True)

ani.save(
    output_video_path,
    writer=animation.FFMpegWriter(
        fps=30,
        codec="libx264",
        extra_args=["-pix_fmt", "yuv420p"],  # compatibilité maximale
    ),
    dpi=150,
)
plt.close(fig)
# %%
