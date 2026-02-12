#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  1 16:36:19 2025

@author: herreman
"""

import numpy as np
import plaque_simu_lib_3zone as psl
from scipy.optimize import fsolve

#géométrie cuve 
Lg=30      #longueur fluide à gauche (en m)
Ld=80        #longueur plaque à droite (en m)
H=10       #hauteur fluide (en m)

#résolution numérique
Mxg=500     #nombre de points dir x fluide partie gauche
Mxd=1200     #nombre de points dir x plaque partie droite
My=100      #nombre de pts dir y vertical

#parametres de l'excitation 
a=1       #amplitude
freq=0.5      #fréquence (en Hz)

#parametres physiques fluide et plaque 
rho=1000    #densité fluide (en kg m^-3)
rhop=850    #densité plaque (en kg m^-3)
h=0.2     #épaisseur plaque (en m)
longueur_flexion = 6  # l_flex= (D/rho g)**(1/4)  (en m)
gravity=9.81   #accélération grav (en m^2 s^{-1})

#calcul de ki
omega=2*np.pi*freq
def disp(k):
    return omega**2-gravity*k*np.tanh(k*H)
      
ki=fsolve(disp,omega**2/gravity)
 

#grouper les paramètres 
param_geo=[Lg,Ld,H]
param_num=[Mxg,Mxd,My]
param_exc=[a,freq]
param_phy=[rho,rhop,h,longueur_flexion,gravity]

#lancer le calcul
maillages,champs = psl.simu_hydro_elastique(param_geo,param_num,param_exc,param_phy)

#déplier les variables
xg,xm,xd,y,X,Y=maillages
phi_mat,iomzeta_g,iomzeta_d,iomzetap,iomkappa=champs    #phi_mat = potentiel, iomzeta=i*omega*zeta (zeta = surface), iomzetap=i*omega*zeta_p (zeta_p  = plaque), iomkappa=i*omega*kappa (kappa = zetap'' courbure)

#%% animation
import matplotlib.pyplot as plt
from matplotlib import animation

fig=plt.figure(figsize=(4,10))
ax = plt.axes()

Ntime=40
omt=np.linspace(0,Ntime-1,Ntime)*2*np.pi/Ntime


def makefig(nn):
    
    n=nn%Ntime
    
    #print(n)

    expfact=np.exp(1j*omt[n])
    
    x_cuve=[-Lg,-Lg,Ld+Lg,Ld+Lg]
    y_cuve=[H,-H,-H,H]
    
    
    x_onde_g=xg
    y_onde_g=np.real(iomzeta_g/(1j*omega)*expfact)
    
    x_plaque=xm
    y_plaque=np.real(iomzetap/(1j*omega)*expfact)
    
    x_onde_d=xd
    y_onde_d=np.real(iomzeta_d/(1j*omega)*expfact)
 
    
    ax.cla()                                         # vider le graphe précédent
    ax.plot(x_cuve,y_cuve,'k:')
    ax.plot(x_onde_g,y_onde_g,'-',label='onde gauche')  
    ax.plot(x_onde_d,y_onde_d,'-',label='onde droite')         
    ax.plot(x_plaque,y_plaque,'-',label='plaque')        
    ax.set_xlabel('x (m)')                              
    ax.set_ylabel('hauteur (m)')
    ax.set_xlim([-Lg-a-0.1,Ld+Lg])
    ax.set_ylim([-1.1*H,1.1*H])
    #ax.set_ylim([-7e-5,7e-5])
    ax.legend(loc='upper right')                      # ajouter la légende 
    #ax.text(0.5,1.15,'t ='+str(np.around(t[n],decimals=1))+' s')     #ajouter un texte 't=...' (on utilise np.around pour arrondir à 1 décimal, sinon le format est moche). 
    ax.grid()


anim=animation.FuncAnimation(fig,makefig,interval=50,frames=range(0,5*Ntime,1),repeat=True)
#anim.save('test'+str(freq)+'Hz_lflex'+str(longueur_flexion)+'.mp4')

#%% onde de pression dynamique sous la plaque, j'anime le membre de droite de l'équation (inertie plaque inclue mais pas important)
# D d^4 zeta / d x^4 + (rho g - rhop h omega^2 ) zeta =  pression dyna (= rho omega^2 phi) 

#maillage xm,pression dyna pdm  (m = milieu)
pdm = rho*omega**2*phi_mat[My,Mxg:(Mxg+Mxd+1)]

plt.figure(num=2,figsize=(10,5))
plt.plot(xm,pdm.real,label='real')
plt.plot(xm,pdm.imag,label='imag')
plt.xlabel('x (m)')
plt.ylabel('pression dyna en Pa, pour onde d\'amplitude 1 m')
plt.legend()
plt.grid()
plt.show()


#%% animation de la pression dynamique

fig=plt.figure(figsize=(10,5))
ax = plt.axes()

Ntime=40
omt=np.linspace(0,Ntime-1,Ntime)*2*np.pi/Ntime


def makefig(nn):
    
    n=nn%Ntime
    
    #print(n)

    expfact=np.exp(1j*omt[n])
    
    pd_instant=np.real(pdm*expfact)
    
    
    ax.cla()                                         # vider le graphe précédent
    ax.plot(xm,pd_instant,label='dynamic pressure')
    ax.set_xlabel('x (m)')                              
    ax.set_ylabel('pression dyna en Pa, pour onde d\'amplitude 1 m')
    ax.set_ylim([-4e4,4e4])
    ax.legend(loc='upper right')                      # ajouter la légende 
    #ax.text(0.5,1.15,'t ='+str(np.around(t[n],decimals=1))+' s')     #ajouter un texte 't=...' (on utilise np.around pour arrondir à 1 décimal, sinon le format est moche). 
    ax.grid()



anim=animation.FuncAnimation(fig,makefig,interval=50,frames=range(0,5*Ntime,1),repeat=True)
