import numpy as np
import matplotlib.pyplot as plt

import baptiste.files.save as sv
import baptiste.display.display_lib as disp


import baptiste.math.fits as fit
import baptiste.math.RDD as rdd

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

params['facq'] = 1000
params['lm'] = 10000
params['fmax_analyse'] = 40


params['L'] = 20
params['theta_f'] = 0

params['n_bins'] = 24
params['nb_blocks'] = 4

params['voie']='Z'
params['num_geo'] = '4-6'
size = data.shape[0]



#%% 4-6
params['voie']='Z'
params['num_geo'] = '4-6'

params['liste_f_garde'] = [[0,22.3],[30.3,32]] #17-14 Z
params['liste_f_2pi'] = [[23,29.5]] #17-14 Z
params['liste_f_moins2pi'] = []#[[33.5,36]] #17-14 Z
params['f_min'] = 0
params['f_max'] = 12
params['select_data'] = True
params['cut_data'] = False

params,f1,histogram1=fct.histo(params,data, add_nom_fig = params['num_geo'])
f1, phase1 = fct.detect(params, f1, histogram1, add_nom_fig = params['num_geo'])
omega1, k1 = fct.omega_k(params, f1, phase1)
save = False
if save :
    sv.save_all_figs(savefolder, params)

#%% 6-12

params['voie']='Z'
params['num_geo'] = '6-12'

params['liste_f_garde'] = [[0,22.3],[30.3,32]]
params['liste_f_2pi'] = [[23,29.5]]
params['liste_f_moins2pi'] = [[33.5,36]]
params['f_min'] = 0
params['f_max'] = 15
params['select_data'] = False
params['cut_data'] = True

params,f2,histogram2=fct.histo(params,data, add_nom_fig = params['num_geo'])
f2, phase2 = fct.detect(params, f2, histogram2, add_nom_fig = params['num_geo'])
omega2, k2 = fct.omega_k(params, f2, phase2 )

phase2_new = phase2 - 4 * np.pi
f2_new = f2

save = False
if save :
    sv.save_all_figs(savefolder, params)

#%% 12-4 FOIREUX 231303 1714

params['voie']='T'
params['num_geo'] = '12-4'

params['liste_f_garde'] = [[0,22.3],[30.3,32]]
params['liste_f_2pi'] = [[23,29.5]]
params['liste_f_moins2pi'] = [[33.5,36]]
params['f_min'] = 9
params['f_max'] = 21
params['select_data'] = False
params['cut_data'] = False

params,f3,histogram3=fct.histo(params,data, add_nom_fig = params['num_geo'])
f3, phase3 = fct.detect(params, f3, histogram3, add_nom_fig = params['num_geo'])
omega3, k3 = fct.omega_k(params, f3, phase3)

# phase3_new = -phase3 + 5 * np.pi
# f3_new = f3

save = False
if save :
    sv.save_all_figs(savefolder, params)
    
#%% test select 6-12

params['liste_f_garde'] = [[0,22.3],[30.3,32]]
params['liste_f_2pi'] = [[23,29.5]]
params['liste_f_moins2pi'] = [[33.5,36]]
params['f_min'] = 0
params['f_max'] = 15
params['select_data'] = False
params['cut_data'] = True

f2_new, phase2_new = fct.select_data(params, f2, phase2)

phase2_new = -phase2_new + 3 * np.pi

plt.figure()
plt.plot(f2_new,phase2_new)

#%% test select 12-4

params['liste_f_garde'] = [[0,22.3],[30.3,32]]
params['liste_f_2pi'] = [[23,29.5]]
params['liste_f_moins2pi'] = [[33.5,36]]
params['f_min'] = 9
params['f_max'] = 21
params['select_data'] = False
params['cut_data'] = True

f3_new, phase3_new = fct.select_data(params, f3, phase3)

phase3_new = -phase3_new + 5 * np.pi

plt.figure()
plt.plot(f3_new,phase3_new)



    
    

#%% FIND ANGLE


# plt.figure()
# plt.plot(f1,phase1/np.pi,label='G4-G6')
# plt.plot(f2_new,phase2_new/np.pi,label='G6-G12')
# # plt.plot(f3,phase3/np.pi,label='G12-G4')
# plt.title('Déphase en fonction de la fréquence')
# plt.xlabel('fréquence (Hz)')
# plt.ylabel(r'$\phi/\pi$')
# plt.legend()

f1_new = f1[:len(f2_new)]
phase1_new = phase1[:len(f2_new)]

phase2_new = phase2_new

plt.figure()
plt.plot(f1_new, phase1_new,label='G4-G6')
plt.plot(f2_new, phase2_new,label='G6-G12')
plt.title('Déphase en fonction de la fréquence')
plt.xlabel('fréquence (Hz)')
plt.ylabel(r'$\phi$')
plt.legend()

theta_f = np.arctan( (1/np.sqrt(3)) * 2 * ( phase2_new/phase1_new + 0.5))

plt.figure()
plt.plot(f1_new, theta_f, 'r+')
plt.title('Angle de propagation en fonction de la fréquence')
plt.xlabel('fréquence (Hz)')
plt.ylabel(r'$\theta$')



params['theta_f'] = theta_f

#%%

omega, k = fct.omega_k(params, f1_new, phase1_new)
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
    popt, pcov = fit.fit(rdd.RDD_pesante_flexion, k, omega, display = True, err = False, nb_param = 2, p0 = [0.3, 1e4], bounds = [0,[1,50e10]], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
    params[str(params['num_fig'][-1])]['hxrho'] = popt[0]
    params[str(params['num_fig'][-1])]['Dsurrho'] = popt[1]
    params[str(params['num_fig'][-1])]['erreurfit'] = pcov
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)
#%%
type_fit = 'power_law_theta_corrige'
if type_fit == 'power_law_theta_corrige':
    params = disp.figurejolie(params = params, nom_fig = 'k_de_omega')
    popt = fit.fit_powerlaw(k[:],omega[:], display = True, xlabel = '', ylabel = '', legend = '')
    params[str(params['num_fig'][-1])]['power'] = popt
    params[str(params['num_fig'][-1])]['data'] = sv.data_to_dict(['k','omega'], [k,omega], data = [k,omega])
    if save :
        sv.save_graph(savefolder, type_fit, params = params)
#%%
save = True
type_fit = 'pesante_flexion_depth_2m'
if True :#type_fit == 'pesante_flexion_depth':
    params = disp.figurejolie(params = params, nom_fig = type_fit)
    popt, pcov = fit.fit(rdd.RDD_pesante_flexion_depth, k, omega, display = True, err = False, nb_param = 2, p0 = [0.3, 1e4], bounds = [0,[0.4,50e10]], zero = True, th_params = False, xlabel = r'k (m$^{-1}$)', ylabel = r'$\omega$')
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
    
      
if save_matlab :
    sv.save_mat(k, savefolder, title = 'k_' + str(params['index']))
    sv.save_mat(k, savefolder, title = 'omega_' + str(params['index']))
    full_params = {}
    for i in params.keys() :
        if not i.isdigit() :
            full_params[i] = params[i]
    sv.save_mat(full_params, savefolder, title = 'parametres_' + str(params['index']))

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
