#%%
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
#%%

def compute_angle(ypx, ysurf_px, dcm_sur_dpx_vertical=47/202, H_cam=32, d_cam_vitre=82):
    """
    L'angle alpha dépend légèrement de la position qu'on regarde 
    dans la cuve (donc dépend de y, en pixels, sur l'image)
    H_cam est en cm
    d_cam_vitre est en cm
    La fonction retourne alpha en radians, pour une hauteur donnée dans l'image
    """    
    return np.arctan(H_cam/(d_cam_vitre+dcm_sur_dpx_vertical*(ysurf_px-ypx)))


def linfitinterp(x,d=dict,plot=False):
    interp_function = d['interp_function']
    tab_y = np.linspace(np.min(d['tab_ymoy_refmanip']), np.max(d['tab_ymoy_refmanip']))
    X,Y = np.meshgrid(x,tab_y)
    Z = interp_function(X,Y)
    if np.sum(np.isnan(Z)==False)==0:
        return np.zeros(2)*np.nan,np.zeros((2,2))*np.nan
    popt,pcov = curve_fit(lambda x,a,b:a*x + b,Y[np.isnan(Z)==False],Z[np.isnan(Z)==False])
    if plot:
        plt.figure()
        plt.plot(Y.flatten(),Z.flatten(),'o')
        plt.plot(Y.flatten(),popt[0]*Y.flatten()+popt[1])
        plt.show()
    return popt,pcov


def compute_aspect_ratio(x,y,d=dict,plot=False):
    popt,_ = linfitinterp(x,d=d,plot=plot)
    dcm_sur_dpx = popt[0]*y+popt[1]
    return dcm_sur_dpx    
# %%
# test :
#print(np.degrees(compute_angle(290, 490))) # haut de la cuve
#print(np.degrees(compute_angle(490, 490))) # bas de la cuve
# correspond bien aux resultats obtenus à la main
# %%
