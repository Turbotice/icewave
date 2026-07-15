import numpy as np
import matplotlib.pyplot as plt
import pickle

import sys
icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

from vasco.tools.clickonfigures import profile_line_on_image_2clicks
from vasco.tools.clickonfigures import get_n_points_onimage
from vasco.tools.clickonfigures import get_n_points
from vasco.tools.clickonfigures import get_n_points_anyfigure



def click_on_fracture_path_plot_time_evol(n_points, uz,
                                          dict_stereo_pivdata, 
                                          fractures_positions_data,
                                          saveplotdata = True,
                                          saveloc = 'C:/Users/Vasco Zanchi/Desktop/presentations/reunions_hebdo/figures_data_article'):
    """
     dans cette partie on veut cliquer sur quelques pixels qui nous intéressent 
     au début, puis on afficher l'évolution temporelle de l'élévation à ce point 
    
    Rq : le tableau idcs_singl_frac est sous la forme [[xind1,yind1],[xind2,yind2],etc.]
    """
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
    if saveplotdata:
        plt.savefig(f'{saveloc}/fig_timeevol_v1.pdf', dpi=600)
        dict_datafig = {}
        dict_datafig['xvals'] = dict_stereo_pivdata['t']
        dict_datafig['matrix_yvals'] = matrix_temp_evol_uz
        dict_datafig['xvals_vlines'] = times_frac_sec_approx
        dict_datafig['idcs_single_frac'] = idcs_single_frac

        pickle.dump(dict_datafig, open(f'{saveloc}/datafig_timeevol_example.pkl', "wb"))
    plt.show(block=True)


    return matrix_temp_evol_uz, times_frac_sec_approx, idcs_single_frac

def clic2extract_amplitude_withref_singlepixel(temp_evol_uz, t_frac_sec_approx, dict_stereo_pivdata, arr_t_frac_sec=None):
    x2plot = dict_stereo_pivdata['t'] - t_frac_sec_approx
    y2plot = temp_evol_uz
    fig, ax = plt.subplots()
    ax.plot(x2plot, y2plot, 'k')
    if type(arr_t_frac_sec)==np.ndarray:
        for i in range(len(arr_t_frac_sec)):
            ax.vlines(arr_t_frac_sec[i] - t_frac_sec_approx, np.min(y2plot), np.max(y2plot), linestyle='--',color='r', alpha=0.2)

    coords = get_n_points_anyfigure(fig=fig, ax=ax, n_points=3, symbol='+')
    x1 = coords[0][0]
    x2 = coords[1][0]
    x3 = coords[2][0]
    y1 = np.interp(x1, x2plot, y2plot)
    y2 = np.interp(x2, x2plot, y2plot)
    y3 = np.interp(x3, x2plot, y2plot)
    
    ya = np.mean([y1, y3])
    ya_std = np.abs(y1 - y3)
    yb = y2

    ymiddle = (ya+yb)/2
    ymiddle_err = ya_std
    
    # maintenant on choisit un seul point pour le max de uz qui nous intéresse, pour en déduire A_c
    fig, ax = plt.subplots()
    ax.plot(x1, y1, '+r')
    ax.plot(x2, y2, '+r')
    ax.plot(x3, y3, '+r')
    ax.plot(x2plot, y2plot, 'k')
    if type(arr_t_frac_sec)==np.ndarray:
        for i in range(len(arr_t_frac_sec)):
            ax.vlines(arr_t_frac_sec[i] - t_frac_sec_approx, np.min(y2plot), np.max(y2plot), linestyle='--',color='r', alpha=0.2)

    coord = get_n_points_anyfigure(fig, ax, n_points=1)
    amplitude = coord[0][1] - ymiddle
    amplitude_err = ymiddle_err

    dict_results = {}
    dict_results['x1'] = x1
    dict_results['x2'] = x2
    dict_results['x3'] = x3
    dict_results['y1'] = y1
    dict_results['y2'] = y2
    dict_results['y3'] = y3
    dict_results['ymiddle'] = ymiddle
    dict_results['ymiddle_err'] = ymiddle_err
    dict_results['amplitude'] = amplitude
    dict_results['amplitude_err'] = amplitude_err

    return dict_results


def oneclick2extract_amplitude_singlepixel(temp_evol_uz, t_frac_sec_approx, dict_stereo_pivdata, arr_t_frac_sec=None):
    x2plot = dict_stereo_pivdata['t']
    y2plot = temp_evol_uz
    fig, ax = plt.subplots(figsize=(18,12))
    ax.plot(x2plot, y2plot, 'k')
    if type(arr_t_frac_sec)==np.ndarray:
        for i in range(len(arr_t_frac_sec)):
            if arr_t_frac_sec[i]==t_frac_sec_approx:
                clr='b'
            else:
                clr='r'
            ax.vlines(arr_t_frac_sec[i], np.min(y2plot), np.max(y2plot), linestyle='--',color=clr, alpha=0.2)
    ax.set_xlim(t_frac_sec_approx-15, t_frac_sec_approx+10)

    coord = get_n_points_anyfigure(fig=fig, ax=ax, n_points=1)
    x_clic = coord[0][0]
    y_clic = coord[0][1]
    y_interp = np.interp(x_clic, x2plot, y2plot)

    amplitude = np.abs(y_interp)

    dict_results = {}
    dict_results['x_clic'] = x_clic
    dict_results['y_clic'] = y_clic
    dict_results['y_interp'] = y_interp
    dict_results['amplitude'] = amplitude

    return dict_results

def extract_amplitude_singlepixel_automatic(temp_evol_uz, t_frac_sec_approx, dict_stereo_pivdata, arr_t_frac_sec=None):
    """
    Pour cette fonction, on a besoin d'un champ detrended, car sinon 
    les valeurs pour l'amplitude sont décalées en fonction du temps
    La valeur qu'on mesure est la valeur de l'amplitude maximale
    pour tout temps < t_frac_sec_approx
    """
    x2plot = dict_stereo_pivdata['t']
    y2plot = temp_evol_uz

    indices_notbroken = np.where(x2plot <= t_frac_sec_approx)
    y2plot_notbroken = y2plot[indices_notbroken]
    x2plot_notbroken = x2plot[indices_notbroken]

    agmx = np.argmax(np.abs(y2plot_notbroken))
    amplitudemax = np.abs(y2plot_notbroken[agmx])
    elevation_amplitudemax = y2plot_notbroken[agmx]
    time_sec_amplitudemax = x2plot_notbroken[agmx]

    fig, ax = plt.subplots(figsize=(18,12))
    ax.plot(x2plot, y2plot, 'k')
    if type(arr_t_frac_sec)==np.ndarray:
        for i in range(len(arr_t_frac_sec)):
            if arr_t_frac_sec[i]==t_frac_sec_approx:
                clr='b'
            else:
                clr='r'
            ax.vlines(arr_t_frac_sec[i], np.min(y2plot), np.max(y2plot), linestyle='--',color=clr, alpha=0.3)
    plt.plot(time_sec_amplitudemax, elevation_amplitudemax, '+r')
    plt.show()

    dict_results = {}
    dict_results['time_amplitudemax'] = time_sec_amplitudemax
    dict_results['indtime_amplitudemax'] = agmx
    dict_results['amplitudemax'] = amplitudemax
    dict_results['elevation_amplitudemax'] = elevation_amplitudemax

    return dict_results


def click2extract_amplitude(matrix_temp_evol_uz, times_frac_sec_approx, dict_stereo_pivdata):
    #%matplotlib qt

    array_amplitudes_frac = np.zeros(len(times_frac_sec_approx))
    array_full_period_sec = np.zeros(len(times_frac_sec_approx))
    array_t1_sec = np.zeros(len(times_frac_sec_approx))
    array_t2_sec = np.zeros(len(times_frac_sec_approx))
    

    for i in range(len(times_frac_sec_approx)):
        x2plot = dict_stereo_pivdata['t'] - times_frac_sec_approx[i]
        y2plot = matrix_temp_evol_uz[i,:]
        coords = get_n_points(x2plot, y2plot, n_points=2)
        y1 = np.interp(coords[0][0], x2plot, y2plot)
        y2 = np.interp(coords[1][0], x2plot, y2plot)

        array_t1_sec[i] = coords[0][0] + times_frac_sec_approx[i]
        array_t2_sec[i] = coords[1][0] + times_frac_sec_approx[i]
        
        delta_t_sec = coords[1][0] - coords[0][0]
        amplitude_picapic = np.abs(y2 - y1)
        
        array_amplitudes_frac[i] = amplitude_picapic/2

        full_period_sec = delta_t_sec*2
        array_full_period_sec[i] = full_period_sec

    return array_amplitudes_frac, full_period_sec, array_t1_sec, array_t2_sec

def threeclick2extract_amplitude(matrix_temp_evol_uz, times_frac_sec_approx, dict_stereo_pivdata):
    #%matplotlib qt

    array_App1 = np.zeros(len(times_frac_sec_approx))
    array_App2 = np.zeros(len(times_frac_sec_approx))
    array_amplitudes_frac_avg = np.zeros(len(times_frac_sec_approx))
    array_amplitudes_frac_err = np.zeros(len(times_frac_sec_approx))
    array_full_period_sec_1 = np.zeros(len(times_frac_sec_approx))
    array_full_period_sec_2 = np.zeros(len(times_frac_sec_approx))
    array_full_period_sec_avg = np.zeros(len(times_frac_sec_approx))
    array_full_period_sec_err = np.zeros(len(times_frac_sec_approx))    
    array_t1_sec = np.zeros(len(times_frac_sec_approx))
    array_t2_sec = np.zeros(len(times_frac_sec_approx))
    array_t3_sec = np.zeros(len(times_frac_sec_approx))
    

    for i in range(len(times_frac_sec_approx)):
        x2plot = dict_stereo_pivdata['t'] - times_frac_sec_approx[i]
        y2plot = matrix_temp_evol_uz[i,:]
        coords = get_n_points(x2plot, y2plot, n_points=3)
        y1 = np.interp(coords[0][0], x2plot, y2plot)
        y2 = np.interp(coords[1][0], x2plot, y2plot)
        y3 = np.interp(coords[2][0], x2plot, y2plot)

        array_t1_sec[i] = coords[0][0] + times_frac_sec_approx[i]
        array_t2_sec[i] = coords[1][0] + times_frac_sec_approx[i]
        array_t3_sec[i] = coords[2][0] + times_frac_sec_approx[i]
        
        
        delta_t_sec_1 = coords[1][0] - coords[0][0]
        delta_t_sec_2 = coords[2][0] - coords[1][0]
        amplitude_picapic1 = np.abs(y2 - y1)
        amplitude_picapic2 = np.abs(y3 - y2)
        array_App1[i] = amplitude_picapic1
        array_App2[i] = amplitude_picapic2
        
        array_amplitudes_frac_avg[i] = np.mean([amplitude_picapic1, amplitude_picapic2])/2
        array_amplitudes_frac_err[i] = np.std([amplitude_picapic2, amplitude_picapic1])/2

        full_period_sec_1 = delta_t_sec_1*2
        full_period_sec_2 = delta_t_sec_2*2

        full_period_sec_avg = np.mean([full_period_sec_1, full_period_sec_2])
        full_period_sec_err = np.std([full_period_sec_1, full_period_sec_2])
        

        array_full_period_sec_1[i] = full_period_sec_1
        array_full_period_sec_2[i] = full_period_sec_2

        array_full_period_sec_avg[i] = full_period_sec_avg
        array_full_period_sec_err[i] = full_period_sec_err
        
    dict_results = {}
    dict_results['array_App1'] = array_App1
    dict_results['array_App2'] = array_App2
    dict_results['array_amplitudes_frac_avg'] = array_amplitudes_frac_avg
    dict_results['array_amplitudes_frac_err'] = array_amplitudes_frac_err
    
    dict_results['array_full_period_sec_1'] = array_full_period_sec_1
    dict_results['array_full_period_sec_2'] = array_full_period_sec_2
    dict_results['array_full_period_sec_avg'] = array_full_period_sec_avg
    dict_results['array_full_period_sec_err'] = array_full_period_sec_err
    dict_results['array_t1_sec'] = array_t1_sec
    dict_results['array_t2_sec'] = array_t2_sec
    dict_results['array_t3_sec'] = array_t3_sec
    

    return dict_results





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