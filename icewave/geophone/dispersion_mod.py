

#%% INITIALISATION

import numpy as np
import matplotlib.pyplot as plt

import baptiste.files.save as sv
import baptiste.display.display_lib as disp


import baptiste.math.fits as fit
import baptiste.math.RDD as rdd
import baptiste.signal_processing.smooth_signal as smooth
import baptiste.tools.tools as tools

import icewave.geophone.fonctions_extract_triangles_flexion as fct

#%% IMPORT FILES

date = '20230313'

# folder ='/run/user/1000/gvfs/smb-share:server=adour.local,share=data/thiou/labshared2/Banquise/Rimouski_2023/Data/Geophones/'+date+'/'
# savefolder = '/run/user/1000/gvfs/smb-share:server=adour.local,share=data/thiou/labshared2/Banquise/Rimouski_2023/Traitements_donnees/baptsarahantonin'

folder = 'W:\Banquise/Rimouski_2023/Data/Geophones\\' +date+'\\'
savefolder = 'W:\Banquise\Rimouski_2023\Traitements_donnees\\baptsarahantonin\\'

special_folder = 'voie_Z_' + date + '_1714\\'
savefolder = savefolder + special_folder
data=fct.import_files(date,3, folder, savefolder)



#%% PARAMETRES

params = {}
params['special_folder'] = ''
params['facq'] = 1000
params['lm'] = 7200
params['fmax_analyse'] = 16


params['L'] = 20
params['theta_f'] = 0

params['n_bins'] = 24
params['nb_blocks'] = 2
params['weight'] = True

params['savefolder'] = savefolder


#%% 4-6
params['voie']='Z'
params['num_geo'] = '4-6'

params['liste_f_garde'] = [[0,26.3], [26.7,30], [79.5, 100]] #[[0,22.2]] #17-14 Z
params['liste_f_2pi'] = [[40,70]] #17-14 Z
params['liste_f_moins2pi'] = [[100,1100]]#[[33.5,36]] #17-14 Z
params['f_min'] = 0
params['f_max'] = 12
params['select_data'] = False
params['cut_data'] = False

params,f1,histogram1=fct.histo(params,data, add_nom_fig = params['num_geo']+ '_' + params['voie'])
f1, phase1 = fct.detect(params, f1, histogram1, add_nom_fig = params['num_geo']+ '_' + params['voie'], save = False)


#%% 6-12

params['voie']='Z'
params['num_geo'] = '6-12'

params['liste_f_garde'] = [[4.7,17]]
params['liste_f_2pi'] = [[0,4.7]]
params['liste_f_moins2pi'] = [[17,36]]
params['f_min'] = 0
params['f_max'] = 90
params['select_data'] = False
params['cut_data'] = False

params,f2,histogram2=fct.histo(params,data, add_nom_fig = params['num_geo']+ '_' + params['voie'])
f2, phase2 = fct.detect(params, f2, histogram2, add_nom_fig = params['num_geo']+ '_' + params['voie'], save = False)


#%% 12-4

params['voie']='Z'
params['num_geo'] = '12-4'

params['liste_f_garde'] = [[0,22.3],[30.3,32]]
params['liste_f_2pi'] = [[23,29.5]]
params['liste_f_moins2pi'] = [[33.5,36]]
params['f_min'] = 0
params['f_max'] = 90
params['select_data'] = False
params['cut_data'] = False

params,f3,histogram3=fct.histo(params,data, add_nom_fig = params['num_geo']+ '_' + params['voie'])
f3, phase3 = fct.detect(params, f3, histogram3, add_nom_fig = params['num_geo']+ '_' + params['voie'], save = True)

   
#%% Select data 4-6

#17-14
# params['liste_f_garde'] = [[79.7,96]]#[[0,26.3], [79.5, 100]] #[[0,22.2]] #17-14 Z
# params['liste_f_2pi'] = [[41,49.3],[53.5,70.8]]#[ [26.7,30]] #17-14 Z
# params['liste_f_moins2pi'] = [[0,29.3],[99.8,104]]#[[140,1100]]#[[33.5,36]] #17-14 Z

#15-13
params['liste_f_garde'] = [[0,2]]
params['liste_f_2pi'] = [[2.9,12.6]]
params['liste_f_moins2pi'] = [[40,50],[99.8,104]]

params['f_min'] = 0
params['f_max'] = 15.8
params['select_data'] = False
params['cut_data'] = False

f1, phase1 = fct.select_data(params, f1, phase1)

phase1 = phase1 + np.pi


plt.figure()
plt.plot(f1, phase1)


#%% Select data 6-12

#17-14
# params['liste_f_garde'] = [[23.7,28.5],[73,79.5]]#[[5,17],[20.5,25.5]]
# params['liste_f_2pi'] = [[0,22.8],[63.5,70.7],[81,91]]#[[40,40]]#[[18.5,20]]
# params['liste_f_moins2pi'] = [[29.7,36]]#[[0,4.7]]

#15-13
params['liste_f_garde'] = [[0.8,2.5]]
params['liste_f_2pi'] = [[40,42.8]]
params['liste_f_moins2pi'] = [[3.9,12.6]]
params['f_min'] = 0
params['f_max'] = 15.8
params['select_data'] = False
params['cut_data'] = False

f2, phase2 = fct.select_data(params, f2, phase2)

phase2 = phase2 - np.pi

plt.figure()
plt.plot(f2,phase2)


#%% Select data 12-4

params['liste_f_garde'] = [[0.8,1.8],[1.9,3.5]]
params['liste_f_2pi'] = [[0,0.5]]
params['liste_f_moins2pi'] = [[4.3,16],[30.3,32]]
params['f_min'] = 1.9
params['f_max'] = 30
params['select_data'] = False
params['cut_data'] = True

f3, phase3 = fct.select_data(params, f3, phase3)

phase3 = phase3 + 2* np.pi

plt.figure()
plt.plot(f3,phase3)
    

#%% FIND ANGLE

save = False

'''Créé un tableau avec les fréquence communes des deux géophones'''

f1_new, f2_new, phase1_new, phase2_new = tools.commun_space(f1, f2, phase1, phase2, pas_f = params['facq']/params['lm'])

# f1_new, f3_new, phase1_new, phase3_new = tools.commun_space(f1, f3, phase1, phase3, pas_f = params['facq']/params['lm'])

#affiche les deux voies (verifier que phase tend vers 0 pour f tend vers 0)
params = disp.figurejolie(params, nom_fig = 'Phase4_6etphase6_12')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\phi$', f1_new, phase1_new, color = 2, exp = False, legend = 'G4-G6')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\phi$', f2_new, phase2_new, color = 3, exp = False, legend = 'G6-G12')
# params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\phi$', f3_new, phase3_new, color = 4, exp = False, legend = 'G12-G4')

#formule pour theta12 en fonction des phases (4-6/6-12)
theta_f = np.arctan( -1 * (1/np.sqrt(3)) * 2 * ( phase2_new/phase1_new + 0.5))

#formule pour theta12 en fonction des phases (4-6/12-4)
# theta_f = np.arctan( (1/np.sqrt(3)) * 2 * ( phase3_new/phase1_new + 0.5))

#formule pour 23 en fonction des phases (6-12/12-4)
# theta_f = np.arctan( (1/np.sqrt(3)) * 2 * ( phase2_new/phase1_new + 0.5))

params = disp.figurejolie(params, nom_fig = 'theta_f_rawdata')
params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\theta$', f1_new, theta_f, color = 5, legend = r'$\theta (f)$ raw')

theta_f = smooth.savgol(theta_f, 5, 2)
# params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\theta$', f1_new, theta_f, color = 2, exp = False, legend = r'$\theta (f)$ smoothed')
# if save :
#     sv.save_graph (params["savefolder"], nom_fig = 'theta_f_rawdata', params = params)
# #on smooth
# theta_f = smooth.savgol(theta_f, 50, 1)

# params = disp.figurejolie(params, nom_fig = 'theta_f_smooth')
# params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'f (Hz)',r'$\theta$', f1_new, theta_f, color = 5)

if save :
    sv.save_graph (params["savefolder"], nom_fig = 'theta_f_smooth', params = params)

params['theta_f'] = theta_f

#On trouve le omega et k
omega, k = fct.omega_k(params, f1_new, phase1_new)

#%% FIT PESANTE INERTIE


# params = disp.figurejolie(params = params, nom_fig = 'omega_k_thetaf_corrige_g4g6')
# params[str(params['num_fig'][-1])]['data'] = disp.joliplot(r'k (m$^{-1}$)', r'$\omega(Hz)$', k, omega, color = 2, exp = False)
# sv.save_graph(savefolder,nom_fig= 'omega_k_thetaf_corrige_g4g6', params = params)

# omega_og, k_og = fct.omega_k(params, f1_new, phase1_new)
# k_og = k_og * np.cos(params['theta_f'])
# params = disp.figurejolie(params = params, nom_fig = 'omega_k_pas_corrige_g4g6')
# params[str(params['num_fig'][-1])]['data'] = disp.joliplot( r'k (m$^{-1}$)', r'$\omega(Hz)$', k_og, omega_og, color = 2, exp = False)
# sv.save_graph(savefolder,nom_fig= 'omega_k_pas_corrige_g4g6', params = params)

save = False
save_matlab = False

type_fit = 'pesante_flexion_thetaf_corrige'


if type_fit == 'pesante_flexion_thetaf_corrige':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    popt, pcov = fit.fit(rdd.RDD_pesante_flexion, k, omega, display = True, err = False, nb_param = 2, p0 = [0.3, 1e4], bounds = [0,[0.4,50e10]], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['hxrho'] = popt[0]
    params[str(params['num_fig'][-1])]['Dsurrho'] = popt[1]
    params[str(params['num_fig'][-1])]['erreurfit'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)
#%%POWER LAW
save = False
type_fit = 'power_law_theta_corrige'
if type_fit == 'power_law_theta_corrige':
    params = disp.figurejolie(params = params, nom_fig = 'k_de_omega')
    popt = fit.fit_powerlaw(k[:],omega[:], display = True, xlabel = '', ylabel = '', legend = '')
    params[str(params['num_fig'][-1])]['power'] = popt
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)
#%%
save = False
type_fit = 'pesante_flexion_depth_2m'
if True :#type_fit == 'pesante_flexion_depth':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    popt, pcov = fit.fit(rdd.RDD_pesante_flexion_depth, k, omega, display = True, err = False, nb_param = 2, p0 = [0.01, 1e4], bounds = [0,[0.02,50e10]], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['hxrho'] = popt[0]
    params[str(params['num_fig'][-1])]['Dsurrho'] = popt[1]
    params[str(params['num_fig'][-1])]['erreurfit'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)
#%%       
type_fit = 'flexion'      
if type_fit == 'flexion':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    # popt, pcov = fit.fit_ransac (k, omega, thresh = 0.1, display = True, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$', newfig = False)
    popt,pcov = fit.fit(rdd.RDD_flexion, k, omega, display = True, err = False, nb_param = 1, p0 = [1e4], bounds = [0,50e10], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['Dsurrho'] = popt
    params[str(params['num_fig'][-1])]['erreurfit'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)
    

save_matlab = False      
if save_matlab :
    sv.save_mat(k, savefolder, title = 'k_' + str(params['index']))
    sv.save_mat(k, savefolder, title = 'omega_' + str(params['index']))
    full_params = {}
    for i in params.keys() :
        if not i.isdigit() :
            full_params[i] = params[i]
    sv.save_mat(full_params, savefolder, title = 'parametres_' + str(params['index']))

#%% depth+fit h
save = False
type_fit = 'flexion+h+H26dm'
if type_fit == 'flexion+h+H26dm':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    popt,pcov = fit.fit(rdd.RDD_hfit_depth, k, omega, display = True, err = False, nb_param = 1, p0 = [0.4], bounds = [0,10], zero = True, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['h'] = popt
    params[str(params['num_fig'][-1])]['err_h'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    plt.grid()
    if save :
        sv.save_graph(savefolder, type_fit, params = params)

h = 0.4
D = h**3 * 6.07e9 / (12 * ( 1 - 0.37**2))
ld = (D / 1000 / 9.81)**0.25
HH = 2.6

plt.axvline(x=1/ld, color = 'k' ) 
plt.annotate(r'$kL_d = 1$', (0.15, 100))


plt.axvline(x=1/HH, color = 'k')
plt.annotate(r'$kH = 1$', (0.41, 100))


#%% Loi de puissance double

popt, pcov = curve_fit(fct, x, y, p0 = p0, bounds= bounds)
x_range = np.linspace(np.min(x), np.max(x), len(x))

disp.joliplot('', '', x, y, color= 13, exp = True, legend = legend_data)
disp.joliplot('xlabel', 'ylabel', x_range, fct(x_range, popt[0], popt[1]), color= 5, exp = False, legend = r'fit : param 1 = ' + str(round(popt[0],4)) + ' param 2 = ' + str(round(popt[1],4)))
if th_params is not False :
    disp.joliplot(xlabel, ylabel, x_range, fct(x_range, th_params[0], th_params[1]), exp = False, legend = r'Theoretical curve', color = 3,zeros = True)
if err :
    plt.fill_between(x_range, fct(x_range, popt[0] + np.sqrt(np.diag(pcov))[0], popt[1] + np.sqrt(np.diag(pcov))[1]),
                     fct(x_range, popt[0] - np.sqrt(np.diag(pcov))[0], popt[1] - np.sqrt(np.diag(pcov))[1]), color = disp.vcolors(2))
    plt.fill_between(x_range, fct(x_range, popt[0] + np.sqrt(np.diag(pcov))[0], popt[1] - np.sqrt(np.diag(pcov))[1]),
                     fct(x_range, popt[0] - np.sqrt(np.diag(pcov))[0], popt[1] + np.sqrt(np.diag(pcov))[1]), color = disp.vcolors(2))



#%%fi capillaire

save = False
type_fit = 'capilaire'
if type_fit == 'capilaire':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    popt,pcov = fit.fit(rdd.RDD_capilaire, k, omega, display = True, err = False, nb_param = 1, p0 = [0.4], bounds = [0,100000], zero = True, th_params = 0.4, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['gamma'] = popt
    params[str(params['num_fig'][-1])]['err_gamma'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)







#%%Fit power law




#%% Fit pour L et T

# select data
params['f_min'] = 2
params['f_max'] = 15
params['select_data'] = False
params['cut_data'] = True

f1_new, phase1_new = fct.select_data(params,f1, phase1)

phase1_new =  phase1_new #+ 4 * np.pi
f2_new, phase2_new = fct.select_data(params,f2, phase2)
phase2_new =  phase2_new #+ 2 * np.pi
f3_new, phase3_new = fct.select_data(params,f3, phase3)
phase3_new =  phase3_new #+ 2 * np.pi
# Angle

plt.figure()
plt.plot(f1_new, phase1_new,label='G4-G6')
plt.plot(f2_new, phase2_new,label='G6-G12')
plt.title('Déphase en fonction de la fréquence')
plt.xlabel('fréquence (Hz)')
plt.ylabel(r'$\phi$')
plt.legend()

plt.figure()
plt.plot(f2_new, phase2_new,label='G6-G12')
plt.plot(f3_new, phase3_new,label='G12-G4')
plt.title('Déphase en fonction de la fréquence')
plt.xlabel('fréquence (Hz)')
plt.ylabel(r'$\phi$')
plt.legend()

theta_f_12 = np.arctan( (1/np.sqrt(3)) * 2 * ( phase2_new/phase1_new + 0.5))

theta_f_23 = np.arctan( (1/np.sqrt(3)) * 2 * ( phase3_new/phase2_new - 0.5))

plt.figure()
plt.plot(f1_new, theta_f_12, 'r+')
plt.title('Angle de propagation 1-2 en fonction de la fréquence')
plt.xlabel('fréquence (Hz)')
plt.ylabel(r'$\theta$')
plt.legend()

plt.figure()
plt.plot(f1_new, theta_f_23, 'r+')
plt.title('Angle de propagation 2-3 en fonction de la fréquence')
plt.xlabel('fréquence (Hz)')
plt.ylabel(r'$\theta$')
plt.legend()


params['theta_f'] = theta_f_12

#E
omega2, k2 = fct.omega_k(params, f2_new, phase2_new)

c2 = omega2/k2

disp.figurejolie()
disp.joliplot('f (Hz)', r'c2 (m.s$^{-1}$)', f2_new, c2, color = 5, params = False)

disp.figurejolie()
disp.joliplot(r'k ($m^{-1}$)', r'$\omega (Hz)$', k2, omega2, color = 4, params = False)


#%% k avec formule antonin
f1_new = f1[:len(f2_new)]
phase1_new = phase1[:len(f2_new)]

params['theta_f'] = 0 


omega1, k1 = fct.omega_k(params, f1_new, phase1_new)
omega2, k2 = fct.omega_k(params, f2_new, phase2_new)
omega = omega1

k = np.sqrt(4/3 * (k1**2 + k2**2 - k1 * k2))

disp.figurejolie()
disp.joliplot(r'k ($m^{-1}$)', r'$\omega (Hz)$', k, omega, color = 4, params = False)

type_fit = 'flexion'      
if type_fit == 'flexion':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    # popt, pcov = fit.fit_ransac (k, omega, thresh = 0.1, display = True, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$', newfig = False)
    popt,pcov = fit.fit(rdd.RDD_flexion, k, omega, display = True, err = False, nb_param = 1, p0 = [1e4], bounds = [0,50e10], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['Dsurrho'] = popt
    params[str(params['num_fig'][-1])]['erreurfit'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)
