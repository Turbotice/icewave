
"""
Created on Fri Jul  3 23:18:59 2026
@author: moreaul
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pickle
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib
get_ipython().run_line_magic('matplotlib', 'qt')
plt.close('all')
# =========================================================
# PARAMÈTRES UTILISATEUR
# =========================================================
year = '2025'
wave_type = 'flexural_wave'
# Choix : ["0211"], ["0212"], ou ["0211","0212"]
dates_to_plot = ["0211", "0212"]
plot_ma = False
window = 2
# =========================================================
# CONFIG LIBRARY
# =========================================================
CONFIG_LIBRARY = {
    ('2025','0212','flexural_wave'): {
        'acquisition_configs': [
            {'acq_num': 0, 'x0': 58, 'tmin': 164.30},
            {'acq_num': 0, 'x0': 85, 'tmin': 307.29},
            {'acq_num': 1, 'x0': 130, 'tmin': 53.68},
            {'acq_num': 1, 'x0': 152, 'tmin': 135.12},
            {'acq_num': 1, 'x0': 170, 'tmin': 235.50},
            {'acq_num': 1, 'x0': 190, 'tmin': 327.4},
            {'acq_num': 1, 'x0': 211, 'tmin': 404.57},
            {'acq_num': 1, 'x0': 231, 'tmin': 492.65},
            {'acq_num': 2, 'x0': 260, 'tmin': 129.66},
            {'acq_num': 2, 'x0': 290, 'tmin': 308.40},
            {'acq_num': 2, 'x0': 327, 'tmin': 414.66},
            {'acq_num': 2, 'x0': 348, 'tmin': 495.54},
            {'acq_num': 3, 'x0': 387, 'tmin': 71.68},
            {'acq_num': 3, 'x0': 410, 'tmin': 146.79},
            {'acq_num': 3, 'x0': 426, 'tmin': 230.82},
            {'acq_num': 3, 'x0': 442, 'tmin': 310.76},
            {'acq_num': 4, 'x0': 460, 'tmin': 185.17},
            {'acq_num': 4, 'x0': 490, 'tmin': 321.03},
            {'acq_num': 4, 'x0': 505, 'tmin': 563.85},
            {'acq_num': 5, 'x0': 525, 'tmin': 46.18},
            {'acq_num': 5, 'x0': 545, 'tmin': 122.85},
            {'acq_num': 5, 'x0': 565, 'tmin': 200.61},
            {'acq_num': 5, 'x0': 585, 'tmin': 322.68},
        ],
        'L': 50
    },
    ('2025','0211','flexural_wave'): {
        'acquisition_configs': [
            {'acq_num': 0, 'x0': 40, 'tmin': 125.7},
            {'acq_num': 0, 'x0': 60, 'tmin': 356.1},
            {'acq_num': 0, 'x0': 90, 'tmin': 494.4},
            {'acq_num': 0, 'x0': 120, 'tmin': 560.5},
            {'acq_num': 1, 'x0': 150, 'tmin': 245.1},
            {'acq_num': 1, 'x0': 170, 'tmin': 329.6},
            {'acq_num': 2, 'x0': 290, 'tmin': 453.4},
            {'acq_num': 3, 'x0': 365, 'tmin': 189.7},
            {'acq_num': 3, 'x0': 385, 'tmin': 303.7},
            {'acq_num': 3, 'x0': 405, 'tmin': 382.4},
            {'acq_num': 3, 'x0': 420, 'tmin': 468.9},
            {'acq_num': 3, 'x0': 445, 'tmin': 549},
            {'acq_num': 4, 'x0': 485, 'tmin': 106},
            {'acq_num': 4, 'x0': 505, 'tmin': 183},
            {'acq_num': 4, 'x0': 525, 'tmin': 265.6},
        ],
        'L': 40
    }
}
# =========================================================
# PATHS
# =========================================================
# Toutes les données (young_modulus*.pkl et disp_curv_*.pkl) sont dans le
# répertoire courant : plus de sous-dossiers {year}_BICWIN/{date}/DAS.
base_path = './'
# =========================================================
# LOAD + COMPUTE
# =========================================================
results = {}
for date in dates_to_plot:
    config = CONFIG_LIBRARY[(year, date, wave_type)]
    acq_list = config['acquisition_configs']
    L = config['L']
    x0_array = np.array([c['x0'] for c in acq_list])
    x_mid = x0_array + L / 2
    path_E = os.path.join(base_path, 'young_modulus.pkl')
    path_E_unc = os.path.join(base_path, 'young_modulus_with_uncertainties.pkl')
    with open(path_E, 'rb') as f:
        young_modulus = pickle.load(f)
    with open(path_E_unc, 'rb') as f:
        young_modulus_unc = pickle.load(f)
    E_raw = np.array(young_modulus['E'])
    E_x = np.array(list(young_modulus['xposition'].values()))
    E_unc_x = np.array(list(young_modulus_unc['xposition'].values()))
    E_unc = np.array(young_modulus_unc['E_uncertainty'])
    idx = np.argsort(E_unc_x)
    E_unc_x = E_unc_x[idx]
    E_unc = E_unc[idx]
    E_func = lambda x: np.interp(x, E_x, E_raw)
    E_unc_func = lambda x: np.interp(x, E_unc_x, E_unc)
    E_cfg = E_func(x_mid)
    E_unc_cfg = E_unc_func(x_mid)
    results[date] = {
        "x_mid": x_mid,
        "E": E_cfg,
        "E_unc": E_unc_cfg
    }
# =========================================================
# PLOT
# =========================================================
plt.figure(figsize=(10,5))
for date in dates_to_plot:
    r = results[date]
    plt.plot(r["x_mid"], r["E"], '-o', label=f'E {date}')
    plt.fill_between(
        r["x_mid"],
        r["E"] - r["E_unc"],
        r["E"] + r["E_unc"],
        alpha=0.2
    )
plt.xlabel("x-position (m)")
plt.ylabel("Young's modulus E")
plt.title(f"E vs x (comparison {year})")
plt.grid()
plt.legend()
plt.tight_layout()
plt.show()
# =========================================================
# MODÈLE DE DISPERSION (Stein / Squire)
# =========================================================
rho_ice = 910
c_w = 1480
rho_w = 1000
nu = 0.3
def wavenumbers_stein_squire(rho_ice, h, H, E, nu, freq, c_w, rho_w, equation):
    g = 9.81
    D = E * h**3 / (12 * (1 - nu**2))
    k = np.linspace(1e-12, 5, 100000)
    idx_zero = np.zeros(len(freq))
    flag = 0
    for kf in range(len(freq)):
        omeg = 2 * np.pi * freq[kf]
        if omeg == 0:
            flag = 1
            idx_flag = kf
        else:
            cph = omeg / k
            if equation == 'stein':
                func = rho_w / D * (g - omeg / np.lib.scimath.sqrt((1 / cph) ** 2 - (1 / c_w) ** 2))
                func -= h * omeg**2 * rho_w / D
                func += (omeg / cph) ** 4
            elif equation == 'squire':
                coth = 1 / np.tanh(k * H)
                func = omeg**2 * (k * h * rho_ice / rho_w + coth) - D * k**5 / rho_w - g * k
            else:
                raise ValueError("equation must be 'stein' or 'squire'")
            func[np.imag(func) != 0] = -1
            func = func.real
            idx_zero[kf] = np.where(np.diff(np.signbit(func)))[0][0]
    idx_zero = idx_zero.astype(int)
    k_QS = k[idx_zero]
    if flag:
        k_QS[idx_flag] = 0
    return k_QS
# =========================================================
# BOUCLE INVERSION PAR DATE
# Incertitude : Option B - propagation analytique
#   sigma_h = sigma_res / sqrt(sum((dk/dh)_i^2))
# =========================================================
# pas de différence finie pour dk/dh (en mètres d'épaisseur)
DELTA_H = 0.005
thickness_results = {}
for date in dates_to_plot:
    config = CONFIG_LIBRARY[(year, date, wave_type)]
    acq_list = config['acquisition_configs']
    r = results[date]
    x_mid_all = r["x_mid"]
    E_all = r["E"]
    E_unc_all = r["E_unc"]
    # --- load dispersion ---
    # Fichier par date, directement dans le répertoire courant
    # (plus de sous-dossiers {year}_BICWIN/{date}/DAS).
    disp_path = os.path.join(base_path, f"disp_curv_{date}.pkl")
    with open(disp_path, "rb") as f:
        disp = pickle.load(f)
    freq_dict = disp["freq"]
    kQS_dict = disp["k_QS"]
    thickness_list = []
    x_list = []
    for i, cfg in enumerate(acq_list):
        x_mid = x_mid_all[i]
        E_mid = E_all[i]
        sigma_E_mid = E_unc_all[i]
        freq = freq_dict[(cfg["acq_num"], cfg["x0"], cfg["tmin"])]
        k_obs = kQS_dict[(cfg["acq_num"], cfg["x0"], cfg["tmin"])]
        h_vals = np.arange(0.20, 1.01, 0.01)
        errors = []
        best_h = None
        best_err = np.inf
        for h in h_vals:
            k_model = wavenumbers_stein_squire(
                rho_ice, h, 1.0, E_mid, nu,
                freq, c_w, rho_w, equation='squire'
            )
            err = np.linalg.norm(k_model - k_obs)
            errors.append(err)
            if err < best_err:
                best_err = err
                best_h = h
        errors = np.array(errors)
        # =====================================================
        # INCERTITUDE - Option B : propagation analytique
        # =====================================================
        N = len(freq)
        # résidus au meilleur h
        k_model_best = wavenumbers_stein_squire(
            rho_ice, best_h, 1.0, E_mid, nu,
            freq, c_w, rho_w, equation='squire'
        )
        residuals = k_model_best - k_obs
        # écart-type des résidus (1 paramètre ajusté -> dof = N-1)
        dof = max(N - 1, 1)
        sigma_res = np.sqrt(np.sum(residuals**2) / dof)
        # sensibilité dk/dh par différence finie centrée
        # (on protège les bornes de h_vals)
        h_plus = min(best_h + DELTA_H, h_vals[-1])
        h_minus = max(best_h - DELTA_H, h_vals[0])
        actual_delta = (h_plus - h_minus) / 2
        k_model_plus = wavenumbers_stein_squire(
            rho_ice, h_plus, 1.0, E_mid, nu,
            freq, c_w, rho_w, equation='squire'
        )
        k_model_minus = wavenumbers_stein_squire(
            rho_ice, h_minus, 1.0, E_mid, nu,
            freq, c_w, rho_w, equation='squire'
        )
        dk_dh = (k_model_plus - k_model_minus) / (2 * actual_delta)
        denom = np.sum(dk_dh**2)
        if denom > 0:
            sigma_h_res = sigma_res / np.sqrt(denom)
        else:
            sigma_h_res = 0.0
        # =====================================================
        # Propagation de l'incertitude sur E vers h
        # E est tenu fixe pendant l'ajustement de h : si E est en réalité
        # entaché d'une erreur sigma_E, le h qui minimise le résidu se
        # déplace en conséquence. Par le théorème des fonctions implicites
        # appliqué à l'équation normale du moindre carré
        #   d/dh [ sum_i (k_model_i(h,E) - k_obs_i)^2 ] = 0 ,
        # on obtient au premier ordre :
        #   dh/dE = - sum_i (dk_i/dh * dk_i/dE) / sum_i (dk_i/dh)^2
        # =====================================================
        DELTA_E_REL = 0.01  # perturbation relative de E pour la différence finie
        delta_E = DELTA_E_REL * E_mid
        k_model_Eplus = wavenumbers_stein_squire(
            rho_ice, best_h, 1.0, E_mid + delta_E, nu,
            freq, c_w, rho_w, equation='squire'
        )
        k_model_Eminus = wavenumbers_stein_squire(
            rho_ice, best_h, 1.0, E_mid - delta_E, nu,
            freq, c_w, rho_w, equation='squire'
        )
        dk_dE = (k_model_Eplus - k_model_Eminus) / (2 * delta_E)
        dh_dE = -np.sum(dk_dh * dk_dE) / denom if denom > 0 else 0.0
        sigma_h_from_E = abs(dh_dE) * sigma_E_mid
        # incertitude totale sur h : résidus de l'ajustement + incertitude sur E
        # (sources supposées indépendantes l'une de l'autre)
        sigma_h = np.sqrt(sigma_h_res**2 + sigma_h_from_E**2)
        # incertitude symétrique (1 sigma), par construction de la méthode
        h_err_low = sigma_h
        h_err_high = sigma_h
        thickness_list.append((best_h, h_err_low, h_err_high, sigma_h_res, dh_dE))
        x_list.append(x_mid)
    thickness_results[date] = {
        "x": np.array(x_list),
        "h": np.array([t[0] for t in thickness_list]),
        "h_low": np.array([t[1] for t in thickness_list]),
        "h_high": np.array([t[2] for t in thickness_list]),
        "sigma_h_res": np.array([t[3] for t in thickness_list]),
        "dh_dE": np.array([t[4] for t in thickness_list]),
    }
# =========================================================
# PLOT FINAL
# =========================================================
plt.figure(figsize=(9,5))
for date in dates_to_plot:
    r = thickness_results[date]
    x = r["x"]
    h = r["h"]
    h_low = r["h_low"]
    h_high = r["h_high"]
    plt.plot(x, h, '-o', label=date)
    plt.fill_between(
        x,
        h - h_low,
        h + h_high,
        alpha=0.2
    )
plt.xlabel("Position (m)")
plt.ylabel("Ice thickness h (m)")
plt.title("Ice thickness inversion with uncertainty (propagation, Option B)")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
# =========================================================
# MODULE DE FLEXION D = E h^3 / (12 (1-nu^2))
# Propagation d'incertitude à partir de sigma_E et sigma_h, EN TENANT
# COMPTE de la corrélation entre E et h : h est ajusté à E fixé, donc une
# erreur sur E se répercute aussi sur h via dh/dE (cf. boucle d'inversion
# ci-dessus). On écrit :
#   D_hat - D_true ~= (dD/dE + dD/dh * dh/dE) * eps_E + dD/dh * eps_h_res
# où eps_E (bruit sur E, variance sigma_E^2) et eps_h_res (bruit résiduel
# propre à l'ajustement de h, variance sigma_h_res^2) sont indépendants.
# Dérivées analytiques : dD/dE = D/E , dD/dh = 3D/h.
# =========================================================
flexural_results = {}
for date in dates_to_plot:
    rE = results[date]
    rh = thickness_results[date]
    E_vals = rE["E"]
    E_unc_vals = rE["E_unc"]
    h_vals = rh["h"]
    sigma_h_res_vals = rh["sigma_h_res"]
    dh_dE_vals = rh["dh_dE"]
    x_vals = rh["x"]  # même grille (x_mid) que rE["x_mid"], par construction
    D_vals = E_vals * h_vals**3 / (12 * (1 - nu**2))
    dDdE = D_vals / E_vals
    dDdh = 3 * D_vals / h_vals
    eff_dDdE = dDdE + dDdh * dh_dE_vals
    var_D = (eff_dDdE**2) * E_unc_vals**2 + (dDdh**2) * sigma_h_res_vals**2
    D_unc_vals = np.sqrt(var_D)
    flexural_results[date] = {
        "x": x_vals,
        "D": D_vals,
        "D_unc": D_unc_vals
    }
# =========================================================
# PLOT D
# =========================================================
plt.figure(figsize=(9,5))
for date in dates_to_plot:
    r = flexural_results[date]
    plt.plot(r["x"], r["D"], '-o', label=date)
    plt.fill_between(
        r["x"],
        r["D"] - r["D_unc"],
        r["D"] + r["D_unc"],
        alpha=0.2
    )
plt.xlabel("Position (m)")
plt.ylabel("Flexural rigidity modulus D (J)")
plt.title("Flexural rigidity modulus")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()