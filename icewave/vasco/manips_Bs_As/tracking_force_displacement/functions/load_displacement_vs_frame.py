#%%
import numpy as np
import matplotlib.pyplot as plt
from utils import *
import load_matdata_JZmethod


def load_disp_vs_frame(path2dailydir='D:/manips_BsAs/Basler_images/1119/serie2',
                        acqname_found='acq2_expo2000_facq50', method='JZ',
                       n_passages=1,i0=2237, refimg=2237, W=256,
                         xpx_vis=687, ypx_vis=482, xpx_vis_ref=303, ypx_vis_ref=482,
                         use_ref=1, frame_force_application=2112, frame_frac=2282,
                         dmm=0.0438116100766703,dpx=1):
        
    if method=='PIV':
        # deux facons d'importer les données en fonction de la structure du fichier .mat
        # la fonction load_displacement_data_PIV prend cela en compte et s'adapte pour 
        # ouvrir le fichier et renvoyer le tableau u
        
        matfile = f'{path2dailydir}/{acqname_found}/PIV_processed_{n_passages}passages_i0{i0}_refimg{refimg}_W{W}.mat'

        u, v, xpix, ypix, concat_with_zeros = load_displacement_data_PIV(matfile)

        xind = (np.abs(xpix - xpx_vis)).argmin()
        yind = (np.abs(ypix - ypx_vis)).argmin()

        xind_ref = (np.abs(xpix - xpx_vis_ref)).argmin()
        yind_ref = (np.abs(ypix - ypx_vis_ref)).argmin()


        nb_neighbors2avg = 0

        # average neighbours and smooth vs time :
        #u_field = gaussian_filter1d(np.mean(u[:,yind-2:yind+2,xind-2:xind+2],axis=(1,2)),5)
        if (nb_neighbors2avg!=0):
            u_field = np.nanmean(u[:,yind-nb_neighbors2avg:yind+nb_neighbors2avg,xind-nb_neighbors2avg:xind+nb_neighbors2avg],axis=(1,2))
            u_field_ref = np.nanmean(u[:,yind_ref-nb_neighbors2avg:yind_ref+nb_neighbors2avg,xind_ref-nb_neighbors2avg:xind_ref+nb_neighbors2avg],axis=(1,2))
        else:
            u_field = u[:,yind,xind]
            u_field_ref = u[:,yind_ref,xind_ref]

        if use_ref:
            u_field_relative = u_field - u_field_ref
        else:
            u_field_relative = u_field

        plt.figure()
        plt.plot(np.arange(u.shape[0])+i0,u_field, label='moving point')
        plt.plot(np.arange(u.shape[0])+i0,u_field_ref, label='ref (point immobile)')
        plt.plot(np.arange(u.shape[0])+i0,u_field_relative, label='relative motion')

        plt.xlim(i0,frame_frac)
        #plt.ylim(np.min(v_field)/4,1.5 * np.max(v_field))
        plt.legend()
        plt.show()


        plt.figure()
        plt.plot(np.arange(u.shape[0])+i0,u_field * dmm/dpx, label='moving point')
        plt.plot(np.arange(u.shape[0])+i0,u_field_ref* dmm/dpx, label='ref (point immobile)')
        plt.plot(np.arange(u.shape[0])+i0,u_field_relative* dmm/dpx, label='relative motion')

        plt.xlim(i0,frame_frac)
        #plt.ylim(np.min(v_field)/4,1.5 * np.max(v_field))
        plt.xlabel('t (frames) ')
        plt.ylabel('displacement (mm)')
        plt.legend()
        plt.show()

        # save pkl for displacement vs frames

        dict_results = {}
        dict_results['method'] = method
        dict_results['u_field'] = u_field
        dict_results['u_field_mm'] = u_field * dmm/dpx
        dict_results['u_field_ref'] = u_field_ref
        dict_results['u_field_ref_mm'] = u_field_ref * dmm/dpx
        dict_results['u_field_relative'] = u_field_relative
        dict_results['u_field_relative_mm'] = u_field_relative * dmm/dpx
        dict_results['i0'] = i0
        dict_results['refimg'] = refimg
        dict_results['frame_frac'] = frame_frac
        dict_results['frame_force_application'] = frame_force_application

        dict_results['xind'] = xind
        dict_results['yind'] = yind
        dict_results['xind_ref'] = xind_ref
        dict_results['yind_ref'] = yind_ref
        dict_results['dmm_sur_dpx'] = dmm/dpx
        dict_results['xpix'] = xpix
        dict_results['ypix'] = ypix
        dict_results['u'] = u
        dict_results['v'] = v
        dict_results['frames_piv'] = np.arange(u.shape[0])+i0


    elif method == 'JZ':
        #matfile = f'{path2dailydir}/{acqname_found}/tracking_methodJZ/vasco_exemple_loc.mat'
        matfile = f'{path2dailydir}/{acqname_found}/relative_displacement_cylinder.mat'
        u_field_relative = load_matdata_JZmethod.load_displacement_px(path2matfile=matfile)
        
        u_field_relative = u_field_relative - u_field_relative[i0]

        plt.figure()
        plt.plot(np.arange(len(u_field_relative)),u_field_relative, label='relative motion')

        plt.xlim(i0,frame_frac)
        #plt.ylim(np.min(v_field)/4,1.5 * np.max(v_field))
        plt.legend()
        plt.show()

        dict_results = {}
        dict_results['method'] = method
        dict_results['u_field_relative'] = u_field_relative
        dict_results['u_field_relative_mm'] = u_field_relative * dmm/dpx
        dict_results['i0'] = i0
        dict_results['frame_frac'] = frame_frac
        dict_results['frame_force_application'] = frame_force_application    
        dict_results['dmm_sur_dpx'] = dmm/dpx
        dict_results['frames_piv'] = np.arange(len(u_field_relative)) # pour l'instant on appelle ça frames_piv, mais juste pour pas faire beuguer les autres fonctions

    return dict_results
