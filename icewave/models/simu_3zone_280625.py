#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 27 11:41:56 2025

@author: herreman
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 18 14:42:42 2025

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
longueur_flexion = 2  # l_flex= (D/rho g)**(1/4)  (en m)
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



#%% tester le calcul des coeffients incidents, ref, transmis, ...
wave_prop,plate_prop=psl.acoefficient(param_geo, param_num, param_exc, param_phy, maillages, champs)
[ai,ar,at,ag_max,ad_max]=wave_prop
[ap_mean,ap_max,ap_elast_mean,ap_elast_max,kappa_max]=plate_prop
print('Wave parameters:')
print('incident amplitude ai :',ai)
print('reflected amplitude ar :',ar)
print('transmitted amplitude after plate at:',at)
print('maximal wave amplitude left of plate aleft_max:',ag_max)
print('maximal wave amplitude left of plate aright_max:',ad_max)
print('\nPlate parameters:')
print('mean plate amplitude:',ap_mean)
print('maximal plate amplitude:',ap_max)
print('mean elastic plate amplitude:',ap_elast_mean)
print('maximal elastic plate amplitude:',ap_elast_max)
print('maximal curvature:',kappa_max)


#%% pour omega fixe et geometry fixe, comment varient les coeffs avec ld la longueur de flexion ?
import time
import datetime

#nombre de ld a varier
num=1000

#grille logarithmique entre 0.02 m et 0.02*(2**15)  m
ld_tab=1e-2*2**np.linspace(1,15,num)

ai_tab=np.zeros(num)
ar_tab=np.zeros(num)
at_tab=np.zeros(num)
ag_max_tab=np.zeros(num)
ad_max_tab=np.zeros(num)
ap_mean_tab=np.zeros(num)
ap_max_tab=np.zeros(num)
ap_elast_mean_tab=np.zeros(num)
ap_elast_max_tab=np.zeros(num)
kappa_max_tab=np.zeros(num)

for i,longueur_flexion in enumerate(ld_tab):
    start = time.time()
    param_phy=[rho,rhop,h,longueur_flexion,gravity]
    #lancer le calcul
    maillages,champs = psl.simu_hydro_elastique(param_geo,param_num,param_exc,param_phy)
    end = time.time()
    print('Elapsed time: ',end - start,'. Calculation ',i,' of ',num-1)
    
    wave_prop,plate_prop=psl.acoefficient(param_geo, param_num, param_exc, param_phy, maillages, champs)
    [ai,ar,at,ag_max,ad_max]=wave_prop
    [ap_mean,ap_max,ap_elast_mean,ap_elast_max,kappa_max]=plate_prop

    ai_tab[i]=ai
    ar_tab[i]=ar
    at_tab[i]=at
    ag_max_tab[i]=ag_max
    ad_max_tab[i]=ad_max
    ap_mean_tab[i]=ap_mean
    ap_max_tab[i]=ap_max
    ap_elast_mean_tab[i]=ap_elast_mean
    ap_elast_max_tab[i]=ap_elast_max
    kappa_max_tab[i]=kappa_max
    
#données/paramètres dans une librairie
data_simu={'param_geo':param_geo,'param_num':param_num,'param_exc':param_num,'param_phy':param_phy,'ld_tab':ld_tab,'ai_tab':ai_tab,'ar_tab':ar_tab,'at_tab':at_tab,'ag_max_tab':ag_max_tab,'ad_max_tab':ad_max_tab,'ap_mean_tab':ap_mean_tab,'ap_max_tab':ap_max_tab,'ap_elast_mean_tab':ap_elast_mean_tab,'ap_elast_max_tab':ap_elast_max_tab,'kappa_max_tab':kappa_max_tab,'ki':ki}    

#enregistrer la librairie sous un nom unique
md=datetime.datetime.utcnow()
name='data_simu_'+md.replace(microsecond=0).isoformat()
np.save(name,data_simu)

#%% visualiser les coefficients incidents/refélchies/transmis (eau après flotteur)

nom='coeffs_irt_(H,Lg,Ld)='+str(H)+','+str(Lg)+','+str(Ld)+'_freq='+str(freq)+'.pdf'

plt.figure(num=1,figsize=(10,5))
plt.semilogx(ld_tab,ai_tab,label='ai')
plt.semilogx(ld_tab,ar_tab,label='ar')
plt.semilogx(ld_tab,at_tab,label='at')
#plt.semilogx(ld_tab,ag_max_tab,label='ag_max')
#plt.semilogx(ld_tab,ad_max_tab,label='ad_max')
plt.xlabel('longueur de flexion (m)')
plt.legend()
plt.grid()
plt.savefig(nom)
plt.show()


#%% visualiser les amplitudes max  dans la plaque et devant et derrière la plaque

nom='ampmax_gpd_(H,Lg,Ld)='+str(H)+','+str(Lg)+','+str(Ld)+'_freq='+str(freq)+'.pdf'

plt.figure(num=2,figsize=(10,5))
plt.semilogx(ld_tab,ag_max_tab,label='ag_max')
plt.semilogx(ld_tab,ad_max_tab,label='ad_max')
plt.semilogx(ld_tab,ap_elast_max_tab,label='ap_elast_max')
plt.semilogx(ld_tab,ap_elast_mean_tab,label='ap_elast_mean')
#plt.semilogx(ld_tab,_tab,label='at')
#plt.semilogx(ld_tab,ag_max_tab,label='ag_max')
plt.xlabel('longueur de flexion (m)')
plt.legend()
plt.grid()
plt.savefig(nom)
plt.show()

#%% même chose mais plutôt les amplides moyennes dans la plaque

nom='ampmean_p_(H,Lg,Ld)='+str(H)+','+str(Lg)+','+str(Ld)+'_freq='+str(freq)+'.pdf'


plt.figure(num=3,figsize=(10,5))
plt.semilogx(ld_tab,ag_max_tab,label='ag_max')
plt.semilogx(ld_tab,ap_mean_tab,label='ap_mean')
#plt.semilogx(ld_tab,ap_max_tab,label='ap_max')
plt.semilogx(ld_tab,ap_elast_mean_tab,label='ap_elast_mean')
#plt.semilogx(ld_tab,_tab,label='at')
#plt.semilogx(ld_tab,ag_max_tab,label='ag_max')
#plt.semilogx(ld_tab,ad_max_tab,label='ad_max')
plt.xlabel('longueur de flexion (m)')
plt.legend()
plt.grid()
plt.savefig(nom)
plt.show()

#%% même graphe, mais en loglog afin d'y voir des lois d'échelles
nom='ampmax_gpd_(H,Lg,Ld)='+str(H)+','+str(Lg)+','+str(Ld)+'_freq='+str(freq)+'loglog.pdf'

droit_x=[0.5,50]
droit_y=[9e-1,9e-2]
droit_xx=[0.5,50]
droit_yy=[9e-1,9e-3]
droit_xxx=[0.5,50]
droit_yyy=[9e-1,(9e-1)*10**(-3/2)]
droit_xxxx=[40,160]
droit_yyyy=[1e-3,1e-3*(4**-4)]
plt.figure(num=4,figsize=(9,12))
plt.loglog(ld_tab,ag_max_tab,label='ag_max')
plt.loglog(ld_tab,ad_max_tab,label='ad_max')
plt.loglog(ld_tab,ap_elast_max_tab,label='ap_elast_max')
plt.loglog(ld_tab,ap_elast_mean_tab,label='ap_elast_mean')
#plt.loglog(droit_x,droit_y,'k:',label='exposant -1/2')
#plt.loglog(droit_xxx,droit_yyy,'k',label='exposant -3/4')
plt.loglog(droit_xx,droit_yy,'k--',label='exposant -1')
plt.loglog(droit_xxxx,droit_yyyy,'k:',label='exposant -4')
#plt.text()
#plt.semilogx(ld_tab,_tab,label='at')
#plt.semilogx(ld_tab,ag_max_tab,label='ag_max')
plt.xlabel('longueur de flexion (m)')
plt.legend()
plt.grid()
plt.savefig(nom)
plt.show()

