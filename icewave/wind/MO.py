"""
Monin - Obukhov stability functions

see book Stull 1988
"""
import numpy as np

def stability(z/L):
    seuil = 1e-2
    stab = ''
    if z/L > 0:
        stab = 'stable'
    elif (z/L<seuil and z/L>-seuil):
        stab = 'neutral'
    else:
        stab = 'unstable'
    return stab

def MO_Businger(z,L):
    """
    Stability function of Businger-Dyer (Stull 9.7.5)

    INPUT:
        z, L
    OUPUT:
        Phi_m, Phi_h
        
    """
        
    Km = 1
    Kh = 1
    stab = stability(z/L)
    if stab=='stable':
        Phi_m = 1 + (4.7*z/L)
        Phi_h = 0.74 + 4.7*z/L # Km/Kh=0.74
    elif stab=='neutral':
        Phi_m = 1
        Phi_h = 0.74 # Km/Kh=0.74
    elif stab=='unstable':
        Phi_m = (1- 15*z/L)**(-1/4)
        Phi_h = 0.74*(1- 9*z/L)**(-1/4) # Km/Kh=0.74

    return Phi_m, Phi_h

def integ_Phi_m(z, z0, L):
    """
    MeanU/u* = integ_phi_m(z,z0,L)
    """
    

    k = 0.41
    stab = stability(z/l)
    if stab=='stable':
        Psi_m = 4.7*z/L
    elif stab=='neutral':
        Psi_m = 0.
    elif stab=='unstable':
        x = (1-(15*z/L))**(1/4)
        Psi_m = -2*np.ln((1+x)/2) - np.ln((1+x**2)/2) + 2*1/(np.tanh(x)) - np.pi/2

    return (1/k)*(np.ln(z/z0) + Psi_m)



