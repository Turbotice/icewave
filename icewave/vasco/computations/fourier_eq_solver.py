#%%
import numpy as np
import matplotlib.pyplot as plt
#%%
# Paramètres physiques
D = 1.203e-6        # diffusivité thermique glace (wiki)
phi_0 = -150.0   # flux au bord (W/m²)
lmbd = 2.25     # lambda glace
T0 = 0      # température imposée à x = h(t)
L = 0.01         # longueur max du domaine (metres)

# h(t) : front mobile
def h(t):
    return L#0.02 + 0.001 * t   # exemple : front se déplace lentement

# Discrétisation
Nx = 100
Nt = 20000
dx = L / Nx
dt = 0.5 * dx**2 / D  # condition de stabilité (schéma explicite)

x = np.linspace(0, L, Nx+1)
t = np.linspace(0, Nt*dt, Nt+1)

# Initialisation
T = np.zeros((Nt+1, Nx+1))
T[0, :] = 0  # condition initiale

# Boucle temporelle
for n in range(0, Nt):
    # Mise à jour des points intérieurs
    for i in range(1, Nx):
        T[n+1, i] = T[n, i] + D * dt / dx**2 * (T[n, i+1] - 2*T[n, i] + T[n, i-1])
    
    # Condition au bord x=0 : dérivée donnée
    T[n+1, 0] = T[n+1, 1] + dx * (phi_0 / lmbd)

    # Condition au bord x=h(t)
    pos_h = h(t[n])
    i_h = int(pos_h / dx)
    if i_h <= Nx:
        T[n+1, i_h] = T0
        T[n+1, i_h+1:] = T0  # optionnel : au-delà du front, on impose T=T0
#%%
%matplotlib qt
# Visualisation
plt.figure(figsize=(8, 5))
for k in range(0, Nt, Nt//10):
    plt.plot(x, T[k, :], label=f"t={t[k]:.2f}s")
plt.xlabel("x (m)")
plt.ylabel("T (°C)")
plt.title("Évolution de la température selon x et t")
plt.legend()
plt.show()

# %%
def compute_statio_temperature(h, phi_0, lambd=2.25, T0=0):
    return T0 + (phi_0/lambd)*h

print(compute_statio_temperature(1e-2,-150))

fluxes = np.linspace(100, 300, 20) * (-1)
thicknesses = np.linspace(1e-2, 0.1, 20)

FF, TT = np.meshgrid(fluxes, thicknesses)
StatioTemperatures = np.zeros_like(FF)

StatioTemperatures = compute_statio_temperature(TT, FF)

plt.figure()
plt.imshow(StatioTemperatures,extent=[np.min(thicknesses), np.max(thicknesses), np.min(fluxes), np.max(fluxes)],aspect='auto')
plt.colorbar()
plt.show()