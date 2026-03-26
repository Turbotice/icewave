#%%
import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# %%
disk = 'D:/'
path2dict = os.path.join(disk, "manips_BsAs", "Summary", "dictionaries_alldata", "all_force_displacement_dicts.pkl")
with open(path2dict, 'rb') as f:
    all_force_displacement_dicts = pickle.load(f)


# %%
def get_epsilondot_from_time_disp_thickness(time_arr, displacement_arr, h_avg, h_std, L):
    popt, pcov = curve_fit(lambda x,a,b:a*x+b, time_arr, displacement_arr)

    slope = popt[0] #slope is the velocity of sollicitation (unit depending on time and displacement units)
    intercept = popt[1]
    err_slope = np.sqrt(np.diag(pcov)[0])
    err_intercept = np.sqrt(np.diag(pcov)[1])

    epsilondot_avg = (6*h_avg*slope)/(L**2)
    epsilondot_err = np.sqrt( ((6*slope/(L**2))**2) *h_std**2 + ((6*h_avg/(L**2))**2) *err_slope**2 )

    print('epsilondot=',epsilondot_avg,'+-',epsilondot_err)
    print('uncertainty of epsilondot due to h :',h_std*6*slope/(L**2))
    print('uncertainty of epsilondot due to fit :',err_slope*6*h_avg/(L**2))

    return epsilondot_avg, epsilondot_err, slope, err_slope, intercept, err_intercept


def compute_epsilonc_epsilondot_atfracture(i, npts2fit=10, all_force_displacement_dicts=all_force_displacement_dicts):
    L = all_force_displacement_dicts[i]['L']
    w = all_force_displacement_dicts[i]['w']
    h_avg = all_force_displacement_dicts[i]['h_avg']
    h_std = all_force_displacement_dicts[i]['h_std']

    idx_frac_common = all_force_displacement_dicts[i]['idx_frac_common']
    #xdata2fit = all_force_displacement_dicts[i]['xdata2fit']

    displacement_meters_comidx = all_force_displacement_dicts[i]['displacement_meters_comidx']

    timearray2fit = np.arange(len(displacement_meters_comidx)) * (1/all_force_displacement_dicts[i]['freq_acq_Hz'])

    disp_max = displacement_meters_comidx[idx_frac_common-1]
    disp_max_err = (displacement_meters_comidx[idx_frac_common-1] - displacement_meters_comidx[idx_frac_common-2])/2
    print('displacement_max =',disp_max,'m')
    # compute the flexural strain just before fracture
    epsilonc_avg = (6*h_avg*disp_max)/(L**2)
    epsilonc_err = np.sqrt( ((6*disp_max/(L**2))**2) *h_std**2 + ((6*h_avg/(L**2))**2) *disp_max_err**2 ) 
    print('epsilonc=',epsilonc_avg,'+-',epsilonc_err)

    print('uncertainty of epsilonc due to h :',h_std*6*disp_max/(L**2))
    print('uncertainty of epsilonc due to disp_max :',disp_max_err*6*h_avg/(L**2))

    ##############################################################
    # Fit pour obtenir epsilondot avant la fracture

    epsilondot_frac_avg, epsilondot_frac_err, slope_frac, err_slope_frac, intercept_frac, err_intercept_frac = get_epsilondot_from_time_disp_thickness(timearray2fit[idx_frac_common-npts2fit:idx_frac_common],
                                                                                         displacement_meters_comidx[idx_frac_common-npts2fit:idx_frac_common],
                                                                                         h_avg,
                                                                                         h_std,
                                                                                         L)
    vz_frac = slope_frac
    vz_frac_err = err_slope_frac

    if not np.isnan(slope_frac):
        plt.figure()
        plt.plot(timearray2fit[:idx_frac_common], displacement_meters_comidx[:idx_frac_common], '.')
        plt.plot(timearray2fit[idx_frac_common - npts2fit : idx_frac_common],timearray2fit[idx_frac_common - npts2fit : idx_frac_common] * slope_frac + intercept_frac, label='fit : displacement speed='+str(np.round(slope_frac*1e3,decimals=3))+'mm/s')
        plt.xlabel('time (s)')
        plt.ylabel('displacement (m)')
        plt.legend()
        plt.show()
    ##############################################################
    # Fit pour obtenir epsilondot sur la même gamme de forces où on a fitté le module d'young
    # Rq : on a juste à utiliser ydata2fit qui est déjà "coupé" aux indices où on a fitté E_eff
    displacement_arr = all_force_displacement_dicts[i]['ydata2fit']
    time_arr = np.arange(len(displacement_arr)) * (1/all_force_displacement_dicts[i]['freq_acq_Hz'])


    epsilondot_avg, epsilondot_err, slope, err_slope, intercept, err_intercept = get_epsilondot_from_time_disp_thickness(time_arr, displacement_arr, h_avg,h_std, L)
    vz = slope
    vz_err = err_slope

    if not np.isnan(slope_frac):
        plt.figure()
        plt.plot(time_arr, displacement_arr, '.')
        plt.plot(time_arr,time_arr * slope + intercept, label='fit : displacement speed='+str(np.round(slope*1e3,decimals=3))+'mm/s')
        plt.xlabel('time (s)')
        plt.ylabel('displacement (m)')
        plt.legend()
        plt.show()



    return epsilonc_avg, epsilonc_err, epsilondot_avg, epsilondot_err, epsilondot_frac_avg, epsilondot_frac_err, vz, vz_err, vz_frac, vz_frac_err

# %%
for i in range(len(all_force_displacement_dicts)):
    print('i=',i)
    epsilonc_avg, epsilonc_err, epsilondot_avg, epsilondot_err,epsilondot_frac_avg, epsilondot_frac_err, vz, vz_err, vz_frac, vz_frac_err = compute_epsilonc_epsilondot_atfracture(i)
    all_force_displacement_dicts[i]['epsilonc_avg'] = epsilonc_avg
    all_force_displacement_dicts[i]['epsilonc_err'] = epsilonc_err
    all_force_displacement_dicts[i]['epsilondot_avg'] = epsilondot_avg
    all_force_displacement_dicts[i]['epsilondot_err'] = epsilondot_err
    
    all_force_displacement_dicts[i]['epsilondot_frac_avg'] = epsilondot_frac_avg
    all_force_displacement_dicts[i]['epsilondot_frac_err'] = epsilondot_frac_err
    all_force_displacement_dicts[i]['vz'] = vz
    all_force_displacement_dicts[i]['vz_err'] = vz_err
    all_force_displacement_dicts[i]['vz_frac'] = vz_frac
    all_force_displacement_dicts[i]['vz_frac_err'] = vz_frac_err
    
# save dicts, now containing epsilonc and epsilondot info
with open(path2dict, 'wb') as f:
    pickle.dump(all_force_displacement_dicts, f)

# %%
