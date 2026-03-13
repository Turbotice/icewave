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

def compute_epsilonc_epsilondot_atfracture(i, npts2fit=10, all_force_displacement_dicts=all_force_displacement_dicts):
    L = all_force_displacement_dicts[i]['L']
    w = all_force_displacement_dicts[i]['w']
    h_avg = all_force_displacement_dicts[i]['h_avg']
    h_std = all_force_displacement_dicts[i]['h_std']

    idx_frac_common = all_force_displacement_dicts[i]['idx_frac_common']
    #xdata2fit = all_force_displacement_dicts[i]['xdata2fit']
    #ydata2fit = all_force_displacement_dicts[i]['ydata2fit']

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


    #slope, intercept = np.polyfit(timearray2fit[idx_frac_common - npts2fit : idx_frac_common], ydata2fit[idx_frac_common - npts2fit : idx_frac_common], 1)

    #try:
    #    p, cov = np.polyfit(timearray2fit[idx_frac_common-npts2fit:idx_frac_common],
    #                    displacement_meters_comidx[idx_frac_common-npts2fit:idx_frac_common],
    #                    1, cov=True)
    #
    #    slope, intercept = p
    #    err_slope     = np.sqrt(cov[0,0])
    #    err_intercept = np.sqrt(cov[1,1])
    #except Exception as e:
    #    print(f"Polyfit failed: {e}")
    #    slope, intercept, err_slope, err_intercept = np.nan, np.nan, np.nan, np.nan

    popt, pcov = curve_fit(lambda x,a,b:a*x+b,
                            timearray2fit[idx_frac_common-npts2fit:idx_frac_common],
                             displacement_meters_comidx[idx_frac_common-npts2fit:idx_frac_common])

    slope = popt[0]
    intercept = popt[1]
    err_slope = np.sqrt(np.diag(pcov)[0])
    err_intercept = np.sqrt(np.diag(pcov)[1])

    epsilondot_avg = (6*h_avg*slope)/(L**2)
    epsilondot_err = np.sqrt( ((6*slope/(L**2))**2) *h_std**2 + ((6*h_avg/(L**2))**2) *err_slope**2 )

    print('epsilondot=',epsilondot_avg,'+-',epsilondot_err)
    print('uncertainty of epsilondot due to h :',h_std*6*slope/(L**2))
    print('uncertainty of epsilondot due to fit :',err_slope*6*h_avg/(L**2))

    if not np.isnan(slope):
        plt.figure()
        plt.plot(timearray2fit[:idx_frac_common], displacement_meters_comidx[:idx_frac_common], '.')
        plt.plot(timearray2fit[idx_frac_common - npts2fit : idx_frac_common],timearray2fit[idx_frac_common - npts2fit : idx_frac_common] * slope + intercept, label='fit : displacement speed='+str(np.round(slope*1e3,decimals=3))+'mm/s')
        plt.xlabel('time (s)')
        plt.ylabel('displacement (m)')
        plt.legend()
        plt.show()

    return epsilonc_avg, epsilonc_err, epsilondot_avg, epsilondot_err

# %%
for i in range(len(all_force_displacement_dicts)):
    print('i=',i)
    epsilonc_avg, epsilonc_err, epsilondot_avg, epsilondot_err = compute_epsilonc_epsilondot_atfracture(i)
    all_force_displacement_dicts[i]['epsilonc_avg'] = epsilonc_avg
    all_force_displacement_dicts[i]['epsilonc_err'] = epsilonc_err
    all_force_displacement_dicts[i]['epsilondot_avg'] = epsilondot_avg
    all_force_displacement_dicts[i]['epsilondot_err'] = epsilondot_err
    
# save dicts, now containing epsilonc and epsilondot info
with open(path2dict, 'wb') as f:
    pickle.dump(all_force_displacement_dicts, f)

# %%
