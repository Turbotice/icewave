#%%
import numpy as np
import matplotlib.pyplot as plt
import pickle
from utils import str_to_float, correspond_samplenum_acqnum, load_csv_params_allacq
import os
from scipy.optimize import curve_fit


#########################################
# on entre les parametres par défauts de la fonction definie en dessous

dates_to_correct_fil = ['1020'] # date  à corriger (force vs pixel) car fil extensible
Slopes_correction_fil = np.array([0.15384155])

disk_example = 'D:'
data_csv = load_csv_params_allacq()

data_csv_converted = np.zeros(data_csv.shape, dtype=object)  # Create an empty array with the same shape
for i in range(data_csv.shape[0]):
    for j in range(data_csv.shape[1]):
        data_csv_converted[i, j] = str_to_float(data_csv[i, j])

data_csv_example = data_csv_converted

##########################################
##########################################

def load_force_displacement_curve(idx=0, disk=disk_example, data_csv=data_csv_example,
                                  w_sample=0.04,
                                  dates_to_correct_fil=dates_to_correct_fil,
                                    Slopes_correction_fil=Slopes_correction_fil,
                                      plot=True, savepkl=True):

    
    use_ref = int(data_csv[idx, 20])

    frame_frac = data_csv[idx, 12]
    frame_force_application = data_csv[idx,11] #début de l'application de la force

    frame_ref_cam1 = data_csv[idx, 16]
    frame_ref_cam2 = data_csv[idx, 17]

    date = str(int(data_csv[idx, 0])).zfill(4)
    print(date)
    acq = int(data_csv[idx, 1])
    serie = data_csv[idx, 2] # si pas de serie, mettre serie = None

    if np.isnan(serie):
        serie = None
    else:
        serie = int(serie)

    if np.isnan(frame_ref_cam2):
        pass
    else:
        frameref_diff = frame_ref_cam2 - frame_ref_cam1
        frame_frac += frameref_diff

    i0 = int(data_csv[idx, 7])
    refimg = int(data_csv[idx, 8])
    W = int(data_csv[idx, 9])
    n_passages = int(data_csv[idx, 10])

    maxforcefit = float(data_csv[idx, 21])
    minforcefit = float(data_csv[idx, 27])
    qualite_fit = data_csv[idx, 22]
    freq_acq_Hz = float(data_csv[idx, 23])
    time2break_sec = float(data_csv[idx, 24])
    method = data_csv[idx, 25]
    imposed_defect = int(data_csv[idx, 26])

    force_disp_track_dir = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/"

    if serie==None:
        with open(force_disp_track_dir + f'data_force_{date}_acq{acq}.pkl', 'rb') as file:
            dict_force = pickle.load(file)
            print('dict force loaded')
    else:
        with open(force_disp_track_dir + f'data_force_{date}_serie{serie}_acq{acq}.pkl', 'rb') as file:
            dict_force = pickle.load(file)
            print('dict force loaded')

    if serie==None:
        with open(force_disp_track_dir + f'data_displacement_{date}_acq{acq}.pkl', 'rb') as file:
            dict_displacement = pickle.load(file)
            print('dict displacement loaded')
    else:
        with open(force_disp_track_dir + f'data_displacement_{date}_serie{serie}_acq{acq}.pkl', 'rb') as file:
            dict_displacement = pickle.load(file)
            print('dict displacement loaded')
            
    # plot force vs displacement
    if np.isnan(frame_ref_cam1):
        common, idx_a, idx_b = np.intersect1d(dict_displacement['frames_piv'], dict_force['frames'], return_indices=True)
    else:
        common, idx_a, idx_b = np.intersect1d(dict_displacement['frames_piv'], dict_force['frames']+frameref_diff, return_indices=True)

    if method=='JZ':
        common_old = common
        common = common_old[common>=i0]
        idx_a = idx_a[common_old>=i0]
        idx_b = idx_b[common_old>=i0]
    
    print("Common values:", common)
    print("Indices in a:", idx_a)
    print("Indices in b:", idx_b)


    u0_relative = dict_displacement['u_field_relative_mm'][idx_a][0]
    frame_frac = dict_displacement['frame_frac']

    idx_frac_common = np.where(common==frame_frac)[0][0]
    if plot:
        plt.figure()
        plt.plot(dict_force['force_grams'][idx_b], dict_displacement['u_field_relative_mm'][idx_a] - u0_relative,'.')
        #plt.xlim(0, ) # limites à bien choisir pour ensuite bien fitter
        plt.ylim(-1e-1,1)
        plt.xlabel('force (grams)')
        plt.ylabel('displacement (mm)')
        plt.savefig(force_disp_track_dir + 'f_grams_vs_displacement_mm.pdf',dpi=600)
        plt.show()

    # convert force in newtons and displacement in meters
    force_newtons_comidx = dict_force['force_grams'][idx_b] * 1e-3 * 9.81
    displacement_meters_comidx = 1e-3 * (dict_displacement['u_field_relative_mm'][idx_a] - u0_relative)

    xdata2fit = force_newtons_comidx[:idx_frac_common]
    ydata2fit = displacement_meters_comidx[:idx_frac_common]

    xdata2fit = np.array(xdata2fit, dtype=float)
    ydata2fit = np.array(ydata2fit, dtype=float)

    print('maxforcefit',maxforcefit,type(maxforcefit))
    print(qualite_fit,type(qualite_fit))

    print(minforcefit, type(minforcefit))
    if (not (np.isnan(maxforcefit)))&(np.isnan(minforcefit)):
        mask = xdata2fit <= maxforcefit
        xdata2fit = xdata2fit[mask]
        ydata2fit = ydata2fit[mask]
    elif (not (np.isnan(minforcefit))) & (np.isnan(maxforcefit)):
        mask = xdata2fit >= minforcefit
        xdata2fit = xdata2fit[mask]
        ydata2fit = ydata2fit[mask]
    elif (not (np.isnan(minforcefit))) & (not(np.isnan(maxforcefit))):
        mask = (xdata2fit >= minforcefit) & (xdata2fit <= maxforcefit)
        xdata2fit = xdata2fit[mask]
        ydata2fit = ydata2fit[mask]
    else:
        pass

    try:
        popt, pcov = curve_fit(lambda x,a,b: a*x+b , xdata2fit, ydata2fit)
    except:
        popt, pcov = np.array([np.nan,np.nan]), np.ones((2,2))*np.nan
    slope, intercept = popt[0], popt[1]
    slope_err, intercept_err = np.diag(pcov)[0], np.diag(pcov)[1]
    if plot:
        plt.figure()
        plt.plot(xdata2fit, ydata2fit,'.')
        #plt.xlim(0, ) # limites à bien choisir pour ensuite bien fitter
        #plt.ylim(-1e-1,1)
        plt.plot(xdata2fit, slope*xdata2fit+intercept)
        plt.xlabel('Fz (N)')
        plt.ylabel('displacement (m)')
        plt.title(f'{date}, serie {serie}, acq {acq}')
        plt.savefig(force_disp_track_dir + f'f_N_vs_displacement_m_{date}_serie{serie}_acq{acq}.pdf',dpi=600)
        plt.show()


    fd, sd = correspond_samplenum_acqnum(date=date, acq=acq, serie=serie)
    print(f'acq{acq} and serie{serie} correspond to sample name : M{fd}{sd}')
    print(fd)
    print(sd)
    
    prefix_filename_thickness = f'M{fd}{sd}'


    thickness_path1 = f"{disk}/manips_BsAs/epaisseurs/{date}/{prefix_filename_thickness}.txt"
    thickness_path2 = f"{disk}/manips_BsAs/epaisseurs/{date}/{prefix_filename_thickness}_camera.txt"
    thickness_path3 = f"{disk}/manips_BsAs/epaisseurs/{date}/{prefix_filename_thickness}_post_frac.txt"
    if os.path.exists(thickness_path1):
        thickness_path = thickness_path1
    elif os.path.exists(thickness_path2):
        thickness_path = thickness_path2
    else:
        thickness_path = thickness_path3

    data_thickness = np.loadtxt(thickness_path,skiprows=1)
    print('thicknesses measured for this sample:',data_thickness)
    h_avg_mm = np.mean(data_thickness)
    h_std_mm = np.std(data_thickness)

    h_avg = 1e-3 * h_avg_mm
    h_std = 1e-3 * h_std_mm


    # compute young's modulus with the formula for elastic beam

    w = w_sample
    L = data_csv[idx, 15]

    E_beam = (1/(4*w))*(L/h_avg)**3 * (1/slope)
    E_err_fit = (slope_err/slope) * E_beam # ca serait l'incertitude liée au fit, qui est faible par rapport à celle liée à l'épaisseur

    E_beam_upp = (1/(4*w))*(L/(h_avg-h_std))**3 * (1/slope)

    E_err_h = E_beam_upp - E_beam

    if E_err_fit>E_err_h:
        E_err = E_err_fit
    else:
        E_err = E_err_h

    print('E = ',E_beam*1e-9,' GPa')
    print('E_upp = ',E_beam_upp*1e-9,' GPa')

    print('E_upp - E = ',(E_beam_upp-E_beam)*1e-9,' GPa')

    print(idx)

    # save results in a dict:

    results_dict = {}
    results_dict['date'] = date
    results_dict['acq'] = acq
    results_dict['serie'] = serie
    results_dict['idx_frac_common'] = idx_frac_common
    results_dict['prefix_filename_thickness'] = prefix_filename_thickness
    results_dict['force_newtons_comidx'] = force_newtons_comidx
    results_dict['displacement_meters_comidx'] = displacement_meters_comidx
    results_dict['maxforcefit'] = maxforcefit
    results_dict['minforcefit'] = minforcefit
    results_dict['qualite_fit'] = qualite_fit
    results_dict['xdata2fit'] = xdata2fit
    results_dict['ydata2fit'] = ydata2fit
    results_dict['w'] = w
    results_dict['L'] = L
    results_dict['E'] = E_beam
    results_dict['E_err'] = E_err
    results_dict['formula E'] = 'Elastic beam theory: E = (1/(4*w))*(L/h)**3 * (1/slope)'
    results_dict['slope'] = slope
    results_dict['slope_err'] = slope_err
    results_dict['slope info'] = 'slope of the linear function: displacement = f(force), so E is proportional to 1/slope'
    results_dict['intercept'] = intercept
    results_dict['intercept_err'] = intercept_err
    results_dict['h_avg'] = h_avg
    results_dict['h_std'] = h_std
    results_dict['freq_acq_Hz'] = freq_acq_Hz
    results_dict['time2break_sec'] = time2break_sec
    results_dict['imposed_defect'] = imposed_defect
    
    # save results dict in a pkl file for combined displacement and force data
    if serie==None:
        dict_path = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/results_dict_force_displacement_{date}_acq{acq}.pkl"
    else:
        dict_path = f"{disk}/manips_BsAs/Summary/tracking_force_displacement/results_dict_force_displacement_{date}_serie{serie}_acq{acq}.pkl"
    
    if savepkl:
        with open(dict_path, 'wb') as file:
            pickle.dump(results_dict, file)
            print('pkl file saved')
    else:
        print('pkl file not saved')

    return (results_dict, h_avg, h_std, E_beam, E_err,prefix_filename_thickness,date)

