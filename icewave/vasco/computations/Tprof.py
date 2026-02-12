#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf

# %% ce code calcule la "similary solution" du probleme de stefan pour une glace qui flotte à la surface de l'eau
# eau est à Tm = 0°C (temperature de solidification glace eau douce)
# on veut le profil de temperature dans l'épaisseur en fonction de temps
# air au dessus est à temperature T0

kappa_ice = 1.203e-6 # diffusivité thermique
lambda_ice = 2.25 # conductivité thermique
cp_ice = 2.04 # specific heat

def compute_Tprofile(z,t,T0=-10,Tm=0,lamb=lambda_ice,kappa=kappa_ice,cp=cp_ice):
    eta = z/(2*(kappa*t))
    T = T0 + (Tm-T0)*(erf(eta)/erf(lamb))
    return T


# %%

array_z = np.linspace(0,0.03,100)

# fix t to start :
t_fix = 3600 
# compute thickness ice :
a_t_fix = 2*lambda_ice*np.sqrt(kappa_ice*t_fix)
print(a_t_fix) # trop grand par rapport à manip grenoble

plt.plot(array_z, compute_Tprofile(array_z, t_fix))