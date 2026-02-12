#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 17:47:25 2025

@author: moreaul
"""

import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
import gpxpy
import pandas as pd
from scipy.fft import fft2, fftshift, fftfreq, fft, ifft
from matplotlib.path import Path
from datetime import datetime
from scipy.linalg import svd
from scipy.interpolate import make_interp_spline
import pickle



import numpy as np
import math
#import mplcursors

import time
from matplotlib.path import Path
get_ipython().run_line_magic('matplotlib', 'qt')
date = '0212'
year = '2025'
path2data = f'/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/{year}_BICWIN/{date}/DAS/'

path2values_h = os.path.join(path2data, 'thicknesses.pkl')
path2values_E = os.path.join(path2data, 'young_modulus.pkl')

with open(path2values_h, 'rb') as f:
    thicknesses = pickle.load(f)
    print("Contents of thicknesses.pkl:")
    print(thicknesses)

with open(path2values_E, 'rb') as f:
    young_modulus = pickle.load(f)
    print("\nContents of young_modulus.pkl:")
    print(young_modulus)


import numpy as np

# --- Prepare h data ---
h_raw = np.array(thicknesses['h'])
h_xpos_raw = np.array(list(thicknesses['xposition'].values()))

# Compute moving average (window = 4)
window = 2
h_ma = np.convolve(h_raw, np.ones(window)/window, mode='valid')

# Adjust x-position for moving average (centered)
h_xpos_ma = h_xpos_raw[(window - 1)//2 : -(window//2)] if len(h_xpos_raw) >= window else []

# --- Prepare E data ---
E_raw = np.array(young_modulus['E'][:-1])  # drop last value
E_xpos_raw = np.array(list(young_modulus['xposition'].values())[:-1])

E_ma = np.convolve(E_raw, np.ones(window)/window, mode='valid')
E_xpos_ma = E_xpos_raw[(window - 1)//2 : -(window//2)] if len(E_xpos_raw) >= window else []



fig, ax1 = plt.subplots(figsize=(10, 5))

# --- Option to plot moving averages ---
plot_ma = False  # change to True if you want moving averages

# --- Plot raw thickness (with faint line + crosses) ---
ax1.plot(
    h_xpos_raw, h_raw,
    '-', color='tab:blue', alpha=1, linewidth=1,
    label='Raw h'
)
ax1.plot(
    h_xpos_raw, h_raw,
    '+', color='tab:blue', alpha=1, markersize=6
)

ax1.set_xlabel('x-position (m)')
ax1.set_ylabel('Thickness h (m)', color='tab:blue')
ax1.tick_params(axis='y', labelcolor='tab:blue')

if plot_ma:
    ax1.plot(
        h_xpos_ma, h_ma,
        'o-', color='tab:blue', label='Moving Avg h'
    )

# --- Plot raw Young’s modulus (with faint line + crosses) ---
ax2 = ax1.twinx()
ax2.plot(
    E_xpos_raw, E_raw,
    '-', color='tab:red', alpha=0.3, linewidth=1,
    label='Raw E'
)
ax2.plot(
    E_xpos_raw, E_raw,
    '+', color='tab:red', alpha=0.6, markersize=6
)

ax2.set_ylabel("Young's modulus E (GPa)", color='tab:red')
ax2.tick_params(axis='y', labelcolor='tab:red')

if plot_ma:
    ax2.plot(
        E_xpos_ma, E_ma,
        's--', color='tab:red', label='Moving Avg E'
    )

# --- Title, layout, grid ---
plt.title("Thickness and Young’s Modulus vs x-position")
fig.tight_layout()
ax1.grid(True)
# %%


# --- Combine legends from both axes ---
lines1, labels1 = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
# ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right')

plt.show()




nu = 0.3

# --- Interpolate E onto h_xpos_raw ---
E_interp = np.interp(h_xpos_raw, E_xpos_raw, E_raw)

# Compute D on h_xpos_raw grid
D_raw = E_interp * (h_raw**3) / (12 * (1 - nu**2))

# --- Plot D vs position ---
plt.figure(figsize=(10, 5))
plt.plot(h_xpos_raw, D_raw, 'o-', color='tab:green', label='Flexural Rigidity D')
plt.xlabel('x-position (m)')
plt.ylabel('D (Pa·m³)')
plt.title('Flexural Rigidity vs x-position')
plt.grid(True)
plt.legend()
plt.show()




