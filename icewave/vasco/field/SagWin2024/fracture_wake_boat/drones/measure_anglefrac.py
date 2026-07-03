#%%
import numpy as np
import matplotlib.pyplot as plt

import sys
icewave_path = 'C:/Users/Vasco Zanchi/Documents/git_turbotice/icewave/icewave/'
sys.path.append(icewave_path)

from scipy.interpolate import LinearNDInterpolator
from scipy.ndimage import gaussian_filter1d
from scipy.ndimage import gaussian_filter

from find_wave_front import find_lines, group_lines_to_dict


#%%
def complete_routine_measure_angle_frac_lines(fractures_positions_data:dict, 
                                              tolerance_px=2, smooth=False,
                                              minsize_frac_px=5,
                                              method='linfit',
                                              plot=True):

    Lines = find_lines(fractures_positions_data['binary'], tolerance_px=tolerance_px)

    dict_lines_frac = group_lines_to_dict(lines=Lines)
    if plot:
        for i in range(len(Lines)):
            plt.plot(dict_lines_frac['line_'+str(i)]['x_ind'], dict_lines_frac['line_'+str(i)]['y_ind'])


    # mesure de l'orientation des points de chaque ligne de fracture

    nb_lines_frac = len(dict_lines_frac.keys())

    Xind_arr = []
    Yind_arr = []
    Alpha_arr = []
    Alpha_linfit_arr = []

    smooth = smooth

    for i in range(nb_lines_frac):
        xind_arr = dict_lines_frac['line_'+str(i)]['x_ind']
        yind_arr = dict_lines_frac['line_'+str(i)]['y_ind']

        if smooth:
            xind_arr = gaussian_filter1d(xind_arr, sigma=2)

        slope_arr = (np.diff(yind_arr, prepend=0))/(np.diff(xind_arr, prepend=0))
        alpha_arr = np.arctan2(np.diff(yind_arr, prepend=0), np.diff(xind_arr, prepend=0))

        dict_lines_frac['line_'+str(i)]['slope_arr'] = slope_arr
        dict_lines_frac['line_'+str(i)]['alpha_arr'] = alpha_arr

        slope_linfit, intercept_linfit = np.polyfit(xind_arr,yind_arr,1)
        slope_linfit_arr = np.ones(len(xind_arr)) * slope_linfit
        alpha_linfit_arr = np.arctan(slope_linfit_arr)
        intercept_linfit_arr = np.ones(len(xind_arr)) * intercept_linfit

        
        if len(xind_arr)>=minsize_frac_px:
            Xind_arr.append(list(xind_arr))
            Yind_arr.append(list(yind_arr))
            Alpha_arr.append(list(alpha_arr))
            Alpha_linfit_arr.append(list(alpha_linfit_arr))




    Xind_arr = np.array([item for sublist in Xind_arr for item in sublist])
    Yind_arr = np.array([item for sublist in Yind_arr for item in sublist])
    Alpha_arr = np.array([item for sublist in Alpha_arr for item in sublist])
    Alpha_linfit_arr = np.array([item for sublist in Alpha_linfit_arr for item in sublist])

    if method=='linfit':
        interp_func_angles = LinearNDInterpolator(np.column_stack([Xind_arr, Yind_arr]), Alpha_linfit_arr)
    elif method=='diff':
        interp_func_angles = LinearNDInterpolator(np.column_stack([Xind_arr, Yind_arr]), Alpha_arr)

    Xinterp, Yinterp = np.meshgrid(np.linspace(np.min(Xind_arr),np.max(Xind_arr)), np.linspace(np.min(Yind_arr),np.max(Yind_arr)))
    Alpha_interp = interp_func_angles(Xinterp, Yinterp)
    if plot:
        plt.figure()
        plt.pcolormesh(Xinterp, Yinterp, Alpha_interp * 180/np.pi, shading='auto', cmap='viridis', vmin=0, vmax=90)
        plt.plot(Xind_arr, Yind_arr, '.')
        plt.colorbar()
        plt.show()

    dict_results = {}
    dict_results['Xind_arr'] = Xind_arr
    dict_results['Yind_arr'] = Yind_arr
    dict_results['Alpha_arr'] = Alpha_arr
    dict_results['Alpha_linfit_arr'] = Alpha_linfit_arr
    dict_results['interp_func_angles'] = interp_func_angles
    dict_results['Xinterp'] = Xinterp
    dict_results['Yinterp'] = Yinterp
    dict_results['Alpha_interp'] = Alpha_interp
    
    return dict_lines_frac, dict_results

def projection_coord_oninclinated_line(x_arr:np.ndarray, y_arr:np.ndarray, slope_line:float, intercept_line:float):
    alpha = np.arctan(slope_line)
    def linear_func(x,a=slope_line,b=intercept_line):
        return a*x+b

    xprim_along_line = (linear_func(x_arr)-linear_func(x_arr[0]))/np.sin(alpha) + ((linear_func(x_arr[0]) - y_arr[0]) - (linear_func(x_arr) - y_arr))*np.sin(alpha)

    return xprim_along_line