#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 25 15:28:41 2025

@author: herreman
"""

# librairies
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from scipy.optimize import fsolve


def simu_hydro_elastique(param_geo,param_num,param_exc,param_phy):

    #déplier les paramètres
    Lg,Ld,H=param_geo
    Mxg,Mxd,My=param_num
    a,freq=param_exc
    rho,rhop,h,longueur_flexion,gravity=param_phy

    #constantes
    omega=2*np.pi*freq   #(rad s^{-1})
    D=(longueur_flexion**4)*(rho*gravity)   #module de flexion

    
    #pas d'espaces
    dxg=Lg/Mxg   #zone gauche selon x
    dxd=Ld/Mxd   #zone droite selon x
    dy=H/My      #selon y
    
    #nombre total de pts     
    Ng=Mxg*(My+1)         #nombre de points zone eau gauche et droite
    Nd=(Mxd+1)*(My+1)     #nombre de points zone milieu
    N=Nd+2*Ng+2*(Mxd+1)     #nombre de points total  = 2(Mxd + 1) en plus de 2*Ng+Nd pour iomzeta=alpha et iomkap=beta)


    #calcul de ki
    def disp(k):
        return omega**2-gravity*k*np.tanh(k*H)
        
    ki=fsolve(disp,omega**2/gravity)
    print('ki finite',ki,'compared to deep estimate',omega**2/gravity)

    #maillages 1D
    xg=np.linspace(-Lg,-dxg,Mxg)      
    xm=np.linspace(0,Ld,Mxd+1) 
    xd=np.linspace(Ld+dxg,Ld+Lg,Mxg) 
    y=np.linspace(-H,0,My+1) 

    #init vecteur g
    g=np.zeros(N,dtype=np.cdouble())

    #nombre d'éléments non-nuls dans la matrice A  
    #            G            D                     B                      H                                  Int
    N_nonnuls= 3*(My+1)  + 3*(My+1) + (3*(Mxg-1)+3*(Mxd+1)+3*(Mxg-1)) + (3*(Mxg-1)+4*(Mxd+1)+3*(Mxg-1)) + 5*(Mxg-1+Mxd+1+Mxg-1)*(My-1) 
    
    #          G  D  eq kappa    cond dyna 
    N_nonnuls+=3 +3 +4*(Mxd-1)+5*(Mxd-1)   #pour l'hydro-élasticité

    # Initialiser les tableaux d'indices de ligne (row), colonne (col) et valeurs (val) commes des zeros.
    row=np.zeros(N_nonnuls)
    col=np.zeros(N_nonnuls)
    val=np.zeros(N_nonnuls,dtype=np.cdouble())
    

    #compte
    compte=0
    telG=0
    telB=0
    telH=0
    telD=0
    telint=0
    for j in range(My+1):  #pour j de 0 a My
    
        #zone fluide sous la surface libre, à gauche
        for iG in range(Mxg):      #pour iG de 0 a Mxg-1
        
            ind=iG*(My+1)+j
            indH=ind+1
            indHH=ind+2
            indB=ind-1
            indBB=ind-2
            indG=ind-(My+1)
            indGG=ind-2*(My+1)
            indD=ind+(My+1)
            indDD=ind+2*(My+1)
        
        
            #eqs de bords, attention aux coins
            #si bord G
           
            if iG==0:           #3*(My+1)
                telG+=3
                row[compte],col[compte],val[compte]=ind,ind,-3/(2*dxg) -1j*ki 
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indD,4/(2*dxg)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indDD,-1/(2*dxg)  
                compte=compte+1 
                
                g[ind]=2*1j*a*omega*np.cosh(ki*(y[j]+H))/np.sinh(ki*H)*np.exp(-1j*ki*xg[0])
                    
            
            #si bord B
            if (j==0)and(iG!=0):  #3*(Mxg-1)
                telB+=3
                row[compte],col[compte],val[compte]=ind,ind,-3/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indH,4/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indHH,-1/(2*dy)  
                compte=compte+1 
                
            #si bord H (partie fluide)
            if (j==My)and(iG!=0): 
                telH+=3
                row[compte],col[compte],val[compte]=ind,ind,3/(2*dy)-omega**2/gravity  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indB,-4/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indBB,1/(2*dy)  
                compte=compte+1 
        
            #si point interieur
            if (iG!=0)and(j!=0)and(j!=My):
                telint+=5
                row[compte],col[compte],val[compte]=ind,ind,-2*(1/dxg**2+1/dy**2)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indB,1/dy**2
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indH,1/dy**2  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indG,1/dxg**2
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indD,1/dxg**2  
                compte=compte+1 
            
        #zone fluide sous la plaque à droite
        for iD in range(Mxd+1):      #pour iD de 0 a Mxd
        
            ind=Ng+iD*(My+1)+j
            indH=ind+1
            indHH=ind+2
            indB=ind-1
            indBB=ind-2
            indG=ind-(My+1)
            indGG=ind-2*(My+1)
            indD=ind+(My+1)
            indDD=ind+2*(My+1)
            indalpha=2*Ng+Nd+iD
            indbeta=2*Ng+Nd+(Mxd+1)+iD
        
            
            #si bord B
            if (j==0):
                telB+=3
                row[compte],col[compte],val[compte]=ind,ind,3/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indH,-4/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indHH,1/(2*dy)  
                compte=compte+1 
        
            #si bord H (partie plaque)
            if (j==My): 
                telH+=4
                row[compte],col[compte],val[compte]=ind,ind,3/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indB,-4/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indBB,1/(2*dy)  
                compte=compte+1
                row[compte],col[compte],val[compte]=ind,indalpha,-1  
                compte=compte+1         
            
        
            #si point interieur
            if (j!=0)and(j!=My): 
                telint+=5
                if iD==0: #points intérieurs speciaux gauche
                    row[compte],col[compte],val[compte]=ind,ind,-2*(1/(dxg*dxd)+1/dy**2)  
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indB,1/dy**2
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indH,1/dy**2  
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indG,2/(dxg*(dxg+dxd))
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indD,2/(dxd*(dxg+dxd)) 
                    compte=compte+1 
                    
                elif iD==Mxd: #points intérieurs speciaux droite
                    row[compte],col[compte],val[compte]=ind,ind,-2*(1/(dxg*dxd)+1/dy**2)  
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indB,1/dy**2
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indH,1/dy**2  
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indG,2/(dxd*(dxg+dxd))
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indD,2/(dxg*(dxg+dxd)) 
                    compte=compte+1 
                
                else:                
                    row[compte],col[compte],val[compte]=ind,ind,-2*(1/dxd**2+1/dy**2)  
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indB,1/dy**2
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indH,1/dy**2  
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indG,1/dxd**2
                    compte=compte+1 
                    row[compte],col[compte],val[compte]=ind,indD,1/dxd**2  
                    compte=compte+1 
                    
        #zone fluide sous la surface libre, à droite
        for iG in range(Mxg):      #pour iG de 0 a Mxg-1
         
            ind=Ng+Nd+iG*(My+1)+j
            indH=ind+1
            indHH=ind+2
            indB=ind-1
            indBB=ind-2
            indG=ind-(My+1)
            indGG=ind-2*(My+1)
            indD=ind+(My+1)
            indDD=ind+2*(My+1)
         
         
            #eqs de bords, attention aux coins
            #si bord D
            if iG==Mxg-1:
                telD+=3
                row[compte],col[compte],val[compte]=ind,ind,3/(2*dxg) +1j*ki 
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indG,-4/(2*dxg)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indGG,1/(2*dxg)  
                compte=compte+1 
                     
             
            #si bord B
            if (j==0)and(iG!=Mxg-1):
                telB+=3
                row[compte],col[compte],val[compte]=ind,ind,-3/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indH,4/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indHH,-1/(2*dy)  
                compte=compte+1 
                 
            #si bord H (partie fluide)
            if (j==My)and(iG!=Mxg-1): 
                telH+=3
                row[compte],col[compte],val[compte]=ind,ind,3/(2*dy)-omega**2/gravity  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indB,-4/(2*dy)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indBB,1/(2*dy)  
                compte=compte+1 
         
            #si point interieur
            if (iG!=Mxg-1)and(j!=0)and(j!=My):
                telint+=5
                row[compte],col[compte],val[compte]=ind,ind,-2*(1/dxg**2+1/dy**2)  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indB,1/dy**2
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indH,1/dy**2  
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indG,1/dxg**2
                compte=compte+1 
                row[compte],col[compte],val[compte]=ind,indD,1/dxg**2  
                compte=compte+1 
    
    print('test telG:',telG,'vs. attendu',3*(My+1))
    print('test telD:',telD,'vs. attendu',3*(My+1))
    print('test telB:',telB,'vs. attendu',(3*(Mxg-1)+3*(Mxd+1)+3*(Mxg-1))) 
    print('test telH:',telH,'vs. attendu',(3*(Mxg-1)+4*(Mxd+1)+3*(Mxg-1)))
    print('test telint:',telint,'vs. attendu', 5*(Mxg-1+Mxd+1+Mxg-1)*(My-1))
    
        
    
    
    
    print('compte:',compte)
    print('teltot:',telG+telB+telH+telD+telint)
    print('doit etre',3*(My+1)  + 3*(My+1) + (3*(Mxg-1)+3*(Mxd+1)+3*(Mxg-1)) + (3*(Mxg-1)+4*(Mxd+1)+3*(Mxg-1)) + 5*(Mxg-1+Mxd+1+Mxg-1)*(My-1))
    
    #equation pour plaque
    for iD in range(Mxd+1):      #pour iD de 0 a Mxd
        j=My
        ind=Ng+iD*(My+1)+j
        indalpha=2*Ng+Nd+iD
        indbeta=2*Ng+Nd+(Mxd+1)+iD
    
        if iD==0:
            row[compte],col[compte],val[compte]=indalpha,indbeta,1   #kappa =0 sur ces lignes
            compte=compte+1 
            row[compte],col[compte],val[compte]=indbeta,indbeta+1,4   #kappa' =0 sur ces lignes
            compte=compte+1 
            row[compte],col[compte],val[compte]=indbeta,indbeta+2,-1   #kappa' =0 sur ces lignes
            compte=compte+1 
        elif iD==Mxd:
            row[compte],col[compte],val[compte]=indalpha,indbeta,1   #kappa =0 sur ces lignes
            compte=compte+1 
            row[compte],col[compte],val[compte]=indbeta,indbeta-1,4   #kappa' =0 sur ces lignes
            compte=compte+1 
            row[compte],col[compte],val[compte]=indbeta,indbeta-2,-1   #kappa' =0 sur ces lignes
            compte=compte+1
        else:
            row[compte],col[compte],val[compte]=indalpha,indalpha-1,1/dxd**2  
            compte=compte+1
            row[compte],col[compte],val[compte]=indalpha,indalpha,-2/dxd**2  
            compte=compte+1
            row[compte],col[compte],val[compte]=indalpha,indalpha+1,1/dxd**2  
            compte=compte+1
            row[compte],col[compte],val[compte]=indalpha,indbeta,-1  
            compte=compte+1
          
            row[compte],col[compte],val[compte]=indbeta,indbeta-1,D/dxd**2  
            compte=compte+1
            row[compte],col[compte],val[compte]=indbeta,indbeta,-D*2/dxd**2  
            compte=compte+1
            row[compte],col[compte],val[compte]=indbeta,indbeta+1,D/dxd**2  
            compte=compte+1
            row[compte],col[compte],val[compte]=indbeta,indalpha,-rhop*h*omega**2+rho*gravity  
            compte=compte+1
            row[compte],col[compte],val[compte]=indbeta,ind,-rho*omega**2  
            compte=compte+1
            
    #Creation de la matrice creuse
    A= coo_matrix((val, (row, col)), shape=(N, N),dtype=np.cdouble())

    # conversion 
    A = A.tocsr()         

    # solution avec spsolve
    vect = spsolve(A,g)

    phi=vect[:(2*Ng+Nd)]   #potentiel hydro
    iomzeta_g=omega**2/gravity*phi[[iG*(My+1)+My for iG in range(Mxg)]]   #i*omega*surface eau gauche
    iomzeta_d=omega**2/gravity*phi[[Ng+Nd+iG*(My+1)+My for iG in range(Mxg)]]   #i*omega*surface eau gauche
    
    iomzetap=vect[(2*Ng+Nd):(2*Ng+Nd+(Mxd+1))]  #i*omega*surface plaque droite
    iomkappa=vect[(2*Ng+Nd+(Mxd+1)):]         #i*omega*kappa, courbure droite

    #maxreal
    print('Maximum de la partie réelle de phi:',np.max(np.abs(np.real(phi))))
    print('Maximum de la partie imagi de phi:',np.max(np.abs(np.imag(phi))))

    print('compte=',compte)
    print('N_nonnuls=',N_nonnuls,'\n')
    
    
    x=np.concatenate((xg,xm,xd))
    X,Y=np.meshgrid(x,y)
    
    maillages=[xg,xm,xd,y,X,Y]
    
    phi_mat=np.zeros((My+1,Mxg+Mxd+1+Mxg),np.cdouble()) 
    for i in range(Mxg+Mxd+1+Mxg):
        for j in range(0,My+1):
            ind = i*(My+1)+j
            phi_mat[j,i]=phi[ind]   
   
    champs=[phi_mat,iomzeta_g,iomzeta_d,iomzetap,iomkappa] 
   
    return maillages,champs   

#%% measurements for post-processing

def acoefficient(param_geo,param_num,param_exc,param_phy,maillages,champs):
    
    xg,xm,xd,y,X,Y=maillages
    phi_mat,iomzeta_g,iomzeta_d,iomzetap,iomkappa=champs

    Lg,Ld,H=param_geo
    Mxg,Mxd,My=param_num
    a,freq=param_exc
    rho,rhop,h,longueur_flexion,gravity=param_phy

    #constantes
    omega=2*np.pi*freq   #(rad s^{-1})
    D=(longueur_flexion**4)*(rho*gravity)   #module de flexion

    #pas d'espaces
    dxg=Lg/Mxg   #zone gauche selon x
    dxd=Ld/Mxd   #zone droite selon x
    dy=H/My      #selon y

    #calcul de ki
    def disp(k):
        return omega**2-gravity*k*np.tanh(k*H)
         
    ki=fsolve(disp,omega**2/gravity)
    

    lami=2*np.pi/ki    
    

    #ai et ar
    M_interp=200
    x_interp=np.linspace(-Lg,-Lg+lami,M_interp)
    zeta_g_interp=np.interp(x_interp,xg,(iomzeta_g/(1j*omega)))

    ai=np.abs(np.mean(zeta_g_interp*np.exp(1j*ki*x_interp))/a)
    ar=np.abs(np.mean(zeta_g_interp*np.exp(-1j*ki*x_interp))/a)
    
  
    #at
    M_interp=200
    x_interp=np.linspace(Lg+Ld-lami,Lg+Ld,M_interp)
    zeta_d_interp=np.interp(x_interp,xd,(iomzeta_d/(1j*omega)))

    at=np.abs(np.mean(zeta_d_interp*np.exp(1j*ki*x_interp))/a)
    #ar=np.abs(np.mean(zeta_g_interp*np.exp(-1j*ki*x_interp))/a)
    
    
    ag_max=np.max(np.abs(iomzeta_g/(1j*omega)))/a
    ad_max=np.max(np.abs(iomzeta_d/(1j*omega)))/a
    
    
    #max de zeta_p
    weight=np.ones(Mxd+1)
    weight[0]=1/2
    weight[-1]=1/2
    zetap=iomzetap/(1j*omega)
    zm=np.mean(zetap)
    thetam=(12/Ld**3)*np.sum(zetap*(xm-Ld/2)*dxd*weight)
    
    zetap_elast=zetap-zm-thetam*(xm-Ld/2)
    

    ap_mean=np.mean(np.abs(zetap))/a    
    ap_max=np.max(np.abs(zetap))/a
    ap_elast_mean=np.mean(np.abs(zetap_elast))/a
    ap_elast_max=np.max(np.abs(zetap_elast))/a
    
    kappa_max=np.max(np.abs(iomkappa/(1j*omega)))/a
    
    

    wave_prop=[ai,ar,at,ag_max,ad_max]
    plate_prop=[ap_mean,ap_max,ap_elast_mean,ap_elast_max,kappa_max]

    return wave_prop,plate_prop
    

