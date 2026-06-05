import numpy as np
import matplotlib.pyplot as plt

import sys
icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

from vasco.tools.clickonfigures import profile_line_on_image_2clicks
from vasco.tools.clickonfigures import get_n_points_onimage
from vasco.tools.clickonfigures import get_n_points



def click_on_fracture_path_plot_time_evol(n_points, uz,
                                          dict_stereo_pivdata, 
                                          fractures_positions_data):
    # dans cette partie on veut cliquer sur quelques pixels qui nous intéressent 
    # au début, puis on afficher l'évolution temporelle de l'élévation à ce point 

    binary_frac_positions = np.where(fractures_positions_data['binary'],fractures_positions_data['binary'],np.nan)
    #%matplotlib qt

    coords = get_n_points_onimage(binary_frac_positions,n_points=n_points, symbol='+')
    idcs_single_frac = np.zeros((len(coords),2), dtype=int)
    for i in range(len(coords)):
        idcs_single_frac[i,0] = np.round(coords[i][0]).astype(int)
        idcs_single_frac[i,1] = np.round(coords[i][1]).astype(int)


    matrix_temp_evol_uz = np.zeros((n_points, uz.shape[2]))
    times_frac_sec_approx = np.zeros(n_points)

    for i in range(n_points):
        xind = idcs_single_frac[i][0]
        yind = idcs_single_frac[i][1]
        matrix_temp_evol_uz[i,:] = uz[yind, xind, :]
        times_frac_sec_approx[i] = fractures_positions_data['tmaxs_frac_kapt'][yind, xind]

    #%matplotlib inline
    plt.figure()
    for i in range(n_points):
        plt.plot(dict_stereo_pivdata['t'], matrix_temp_evol_uz[i],'k',alpha=0.2)
        plt.vlines(times_frac_sec_approx[i], ymin=(np.min(matrix_temp_evol_uz)), 
                ymax=np.max(matrix_temp_evol_uz), colors='r',linestyle='--',
                    alpha=0.3)
    plt.show(block=True)

    return matrix_temp_evol_uz, times_frac_sec_approx, idcs_single_frac

def click2extract_amplitude(matrix_temp_evol_uz, times_frac_sec_approx, dict_stereo_pivdata):
    #%matplotlib qt

    array_amplitudes_frac = np.zeros(len(times_frac_sec_approx))

    for i in range(len(times_frac_sec_approx)):
        x2plot = dict_stereo_pivdata['t'] - times_frac_sec_approx[i]
        y2plot = matrix_temp_evol_uz[i,:]
        coords = get_n_points(x2plot, y2plot, n_points=2)
        y1 = np.interp(coords[0][0], x2plot, y2plot)
        y2 = np.interp(coords[1][0], x2plot, y2plot)
        
        delta_t_sec = coords[1][0] - coords[0][0]
        amplitude_picapic = np.abs(y2 - y1)
        
        array_amplitudes_frac[i] = amplitude_picapic/2

        full_period_sec = delta_t_sec

    return array_amplitudes_frac, full_period_sec






def plot_elevation_refnotbroken_and_broken(D, uz, idcs_single_frac, times_frac_sec_approx,
                                           dict_stereo_pivdata, matrix_temp_evol_uz,
                                           fractures_positions_data):


    xdata2plot = idcs_single_frac[:,0]
    ydata2plot = idcs_single_frac[:,1]
    plt.figure()
    plt.plot(xdata2plot, ydata2plot)
    plt.plot(xdata2plot[0], ydata2plot[0], 'or')
    
    slope, intercept = np.polyfit(xdata2plot, ydata2plot, 1)



    alpha = np.arctan(slope)
    deltaD = (slope*xdata2plot[0]+intercept - ydata2plot[0])*np.sin(alpha)


    x_ind_ref_noncasse = xdata2plot[0] - ((D+deltaD)/np.sqrt(1+slope**2))
    y_ind_ref_noncasse = slope*x_ind_ref_noncasse + intercept

    x_ind_ref_noncasse_int = np.round(x_ind_ref_noncasse).astype(int)
    y_ind_ref_noncasse_int = np.round(y_ind_ref_noncasse).astype(int)

    plt.plot(x_ind_ref_noncasse, y_ind_ref_noncasse, 'o')
    plt.plot(np.linspace(xdata2plot[0]-5,xdata2plot[-1]), np.linspace(xdata2plot[0]-5,xdata2plot[-1])*slope+intercept)
    plt.show()

    # maintenant on prend le temps pour lequel il est le plus probable que la fracture s'est passée (en gros) à ce moment
    occurences,binvals = np.histogram(times_frac_sec_approx)
    most_probable_binval = binvals[np.argmax(occurences)]
    time_frac_approx = most_probable_binval

    # on peut donc tracer, pour une valeur de D, l'évolution temporelle de uz associée

    uz_ref_noncassee = uz[y_ind_ref_noncasse_int, x_ind_ref_noncasse_int]
    #%matplotlib qt
    fig, axs = plt.subplots(1,2)
    axs[1].plot(dict_stereo_pivdata['t'], uz_ref_noncassee, label='glace non fracturée')
    for i in range(matrix_temp_evol_uz.shape[0]):
        plt.plot(dict_stereo_pivdata['t'], matrix_temp_evol_uz[i,:],'k',alpha=0.2)

    axs[1].vlines(time_frac_approx, np.min(matrix_temp_evol_uz), np.max(matrix_temp_evol_uz), linestyle='--', color='r', label='instant (approx) de fracture détecté')
    axs[1].plot([],[],'k',alpha=0.2,label='Points situés sur la ligne de fracture')
    axs[1].legend()
    axs[1].set_ylabel('Elevation $u_z$(t) [m]', fontsize=15)
    axs[1].set_xlabel('t [s]', fontsize=15)
    axs[1].set_title('Elevation vs time for points on a given crack')


    axs[0].set_title('Elevation field $u_z$ at time of a given fracture')
    axs[0].imshow(uz[:,:,np.where(dict_stereo_pivdata['t']>=time_frac_approx)[0][0]])
    axs[0].imshow(np.where(fractures_positions_data['binary'],fractures_positions_data['binary'],np.nan),cmap='gray')
    axs[0].plot(x_ind_ref_noncasse_int, y_ind_ref_noncasse_int, 'sr')
    axs[0].plot(idcs_single_frac[:,0], idcs_single_frac[:,1],'r+')

    plt.show()

    return uz_ref_noncassee, x_ind_ref_noncasse_int, y_ind_ref_noncasse_int, time_frac_approx, (x_ind_ref_noncasse,y_ind_ref_noncasse,xdata2plot,ydata2plot,slope,intercept, deltaD)