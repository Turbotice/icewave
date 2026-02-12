# -*- coding: utf-8 -*-
"""
Created on Thu Apr 10 10:09:51 2025

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation

import h5py
import pickle
import os
import glob

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT

# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

#%% Function section

def animation_profile(fig,ax,data,x,phi,nb_frames,time_interval):
    
    """ Create and return animation of a line plot using matplotlib.animation
    Inputs : - fig: matplotlib figure
             - ax: axis object
             - data: numpy array, data to plot #dim0 : space, #dim1 : time
             - x: numpy array, x-axis (space)
             - t: numpy array, (time)
             - nb_frames : int, number of frames to show
             - time_interval : time between two consecutive frames
             
    Output : - ani: matplotlib animation object"""
    
    
    j0 = 0 # initial time index 
    line = ax.plot(x,data[:,j0])[0]
    
    ax.set_title(r'$phi =' + '{:.2f}'.format(phi[j0]) + r' \; \mathrm{(rad)}$')
    ax.set_xlabel(r'$x \; \mathrm{(mm)}$')
    ax.set_ylabel(r'$\xi \; \mathrm{(mm)}$')
    
    # small function to update the current figure 
    def update_profile_plot(frame):
        line.set_xdata(x)
        line.set_ydata(data[:,frame])
        ax.set_title(r'$phi =' + '{:.2f}'.format(phi[frame]) + r' \; \mathrm{(rad)}$')
        return line
    
    # create an animation 
    ani = animation.FuncAnimation(fig=fig, func=update_profile_plot, frames=nb_frames, interval=time_interval)
    plt.show()
    print('Animation computed')

    return ani

def make_movie(V,t,x,y,colormap,cmax,time_interval = 1e3):
    
    """ Create animation of real demodulated fields for successive frequencies 
    Inputs : - velocity field, numpy array, [nx,ny,nf] Time Fourier transform of a 2D field
             - t, numpy array, 1D, array of frequencies
             - colormap, cmap object, colormap used to plot real fields
             - cmax, float, value used to scale colorbar from -cmax to +cmax
             - x, numpy array, array of x coordinates 
             - y, numpy array, array of y coordinates
             
    Outputs : - ani, matplotlib animation which can be saved in a .mp4 format
        """
    nb_frames = len(t) # total number of frames of animation
    fig, ax = plt.subplots(figsize = (12,9))
    
    current_time = t[0]
    print(f't = {current_time}')
    
    field = V[:,:,0]
    # show initial matrix 
    c = ax.imshow(field.T,cmap = colormap,aspect = 'equal', origin = 'lower', interpolation = 'gaussian',
              vmin = -cmax,vmax = cmax,extent = (x.min(),x.max(),y.min(),y.max()))
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$V_x (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)
    
    ax.set_xlabel(r'$x \; \mathrm{(m)}$', labelpad = 5)
    ax.set_ylabel(r'$y \; \mathrm{(m)}$', labelpad = 5)
    
    def update(frame) : 
        current_time = t[frame]
        print(f't = {current_time}')
        
        field = V[:,:,frame]
        c.set_data(field.T)
        ax.set_title('$t = '+ f'{current_time:.3f}'+ '\; \mathrm{(s)}$')
        
        return c
    
    # create an animation 
    ani = animation.FuncAnimation(fig=fig, func=update, frames=nb_frames, interval=time_interval)
    # plt.show()
    print('Animation computed')
    return ani

def animation_complex_field(V,phi,x,y,colormap,cmax,time_interval = 1e3):
    
    nb_frames = len(phi) # total number of frames of animation
    fig, ax = plt.subplots(figsize = (12,9))
    
    current_phi = phi[0]
    print(f'phi = {current_phi}')
    
    field = V*np.exp(-1j*phi[0])
    real_field = np.real(field)
    # show initial matrix 
    c = ax.imshow(real_field.T,cmap = colormap,aspect = 'equal', origin = 'lower', interpolation = 'gaussian',
              vmin = -cmax,vmax = cmax,extent = (x.min(),x.max(),y.min(),y.max()))
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$V_x (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)
    
    ax.set_xlabel(r'$x \; \mathrm{(m)}$', labelpad = 5)
    ax.set_ylabel(r'$y \; \mathrm{(m)}$', labelpad = 5)
    
    def update(frame) : 
        current_phi = phi[frame]
        print(f'phi = {current_phi}')
        
        field = V*np.exp(-1j*phi[frame])
        real_field = np.real(field)
        c.set_data(real_field.T)
        ax.set_title('$phi = '+ f'{current_phi:.2f}'+ '\; \mathrm{(rad)}$')
        
        return c
    
    # create an animation 
    ani = animation.FuncAnimation(fig=fig, func=update, frames=nb_frames, interval=time_interval)
    # plt.show()
    print('Animation computed')
    return ani



def lorentzian(x,x0,alpha):
    y = 1/np.sqrt(1 + ((x - x0)/alpha)**2)
    return y


def set_imshow_extents(c,x,y):
    """ Correct imshow extents so that values are observed in middle of each pixels 
    Inputs : - c, axesImage object (returned by ax.imshow)
             - x, numpy array, array used for setting left/right extents 
             - y, numpy array, array used for setting bottom/top extents"""
      
    shift_x = np.diff(x)[0]*0.5
    shift_y = np.diff(y)[0]*0.5
    extents = (x.min() + shift_x, x.max() + shift_x,
               y.min() + shift_y, y.max() + shift_y)
    c.set_extent(extents)
    
#---------------------------------------------------------------------------------------------------------------
def show_velocity_field(V,x,y,colormap,figname):

    field = V.T
    
    set_graphs.set_matplotlib_param('single')
    fig, ax = plt.subplots(figsize = (12,9))
    c = ax.imshow(field, cmap = parula_map , aspect = 'equal', norm = 'linear', origin = 'lower',interpolation = 'gaussian',
                  extent = (data['x'].min(),data['x'].max(),data['y'].min(),data['y'].max()))
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$V_x (x,y) \; \mathrm{(u.a.)}$',labelpad = 5)
    ax.set_xlabel(r'$x \; \mathrm{(m)}$', labelpad = 5)
    ax.set_ylabel(r'$y \; \mathrm{(m)}$', labelpad = 5)
    
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')

    
#####################################################################
#%% ################### Perform loop ################################
#####################################################################

h = 7.5 # frasil thickness 
date = '2024_07_11'
path2data = f'U:/Aurore_frasil/{date}_e_{h}mm_laser/matData/'

folderlist = glob.glob(f'{path2data}*')

for idx_folder in range(len(folderlist)):
    file2load = glob.glob(f'{folderlist[idx_folder]}/*scaled.mat')[0]
    
    with h5py.File(file2load,'r') as fmat:
        print('Top-level keys : ', list(fmat.keys()))
        data = mat2py.mat_to_dict(fmat['m'],fmat['m'])
    
    print(f'file loaded : {file2load}')
    data = mat2py.transpose_PIVmat_fields(data)
    
    # Define fig_folder and collect variables
    current_path = f'{folderlist[idx_folder]}/'
    fig_folder = f'{current_path}Plots/'
    if not os.path.isdir(fig_folder):
        os.mkdir(fig_folder)
    
    f_ex = data['EXP']['f_ex']
    ID = data['ID']
    ID_txt = f'{h}mm_{ID}'.replace('.','p')
    
    # Plot histogram of velocity 
    
    figname = f'{fig_folder}Histogram_PIV_close_wavemaker_{ID_txt}'
    FT.histogram_PIV(data['Vx']/data['SCALE']['scale_V'],data['PIV_param']['w'],figname)
    
    # Show velocity field 
    frame = 0
    V = data['Vx'][:,:,frame]
    figname = f'{fig_folder}Vx_frame{frame}_{ID_txt}'
    show_velocity_field(V,data['x'],data['y'],parula_map,figname)
    
    # Perform time FFT 
    
    TF_spectrum,freq,FFT_t = FT.temporal_FFT(data['Vx'],data['SCALE']['facq_t'],padding_bool = 1,add_pow2 = 1,output_FFT = True)
    
    fig, ax = plt.subplots()
    set_graphs.set_matplotlib_param('single')
    ax.plot(freq,TF_spectrum)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([1e-5,1e-2])
    
    ax.set_xlabel(r'$f \; \mathrm{(Hz)}$',labelpad = 15)
    ax.set_ylabel(r'$\langle |\hat{V}_x| \rangle _{x,y}(f) \; \mathrm{(u.a.)}$',labelpad = 5)
    
    figname = f'{fig_folder}TF_spectrum_{ID_txt}'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.svg', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    
    # Show demodulated field 
    
    # find peak
    p = np.argmax(TF_spectrum)
    f_max = freq[p]
    print(f'Detected f_ex = {f_max:.2f}')
    
    demod_field = np.mean(data['Vx']*np.exp(1j*2*np.pi*f_max*data['t']),axis = -1)
    real_field = np.real(demod_field)
    
    figname = f'{fig_folder}Demodulated_Vx_{ID_txt}'
    show_velocity_field(real_field, data['x'], data['y'], parula_map, figname)
    
    
    # Perform FFT 2D of demodulated field
    add_pow2 = [1,1]
    facq = (1/data['SCALE']['fx'],1/data['SCALE']['fx']) # acquisition frequency for each direction x and y 
    shift,ky,kx = FT.fft_2D(demod_field,facq,add_pow2)
    
    # find maximum of the 2D FFT
    idx_max = np.argmax(abs(shift).flatten())
    unravel_coords = np.unravel_index(idx_max,shift.shape)
    
    fig, ax = plt.subplots(figsize = (12,9))
    c = ax.imshow(abs(shift).T,cmap = parula_map,aspect = 'equal',norm = 'linear',origin = 'lower',
                  extent = (kx.min(),kx.max(),ky.min(),ky.max()))
    
    ax.plot(kx[unravel_coords[0]],ky[unravel_coords[1]],'ro')
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$\hat{V_x} (x,y) \; \mathrm{(u.a.)}$',labelpad = 1)
    ax.set_xlabel(r'$k_x \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    ax.set_ylabel(r'$k_y \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    
    ax.set_xlim([-100,100])
    ax.set_ylim([-100,100])
    
    figname = f'{fig_folder}Spatial_FFT_{ID_txt}'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    
    # Get profile along ky 
    
    cut = abs(shift)[:,unravel_coords[1]]
    y_exp = (cut - cut.min())/(cut.max() - cut.min())
    
    bounds_kx0 = (kx[unravel_coords[0] - 10],kx[unravel_coords[0] + 10])
    bounds_alpha = (0,1e2)
    bounds_curvefit = ([bounds_kx0[0],bounds_alpha[0]],[bounds_kx0[1],bounds_alpha[1]])
    # fit by a lorentzian
    popt,pcov = scipy.optimize.curve_fit(lambda x,x0,sigma : lorentzian(x, x0, sigma),kx,y_exp,bounds = bounds_curvefit)
    err_coeff = np.sqrt(np.diag(pcov))
    label_fit = r'$\alpha = ' + f'{popt[1]:.1f}' + '\; \mathrm{(m^{-1})}$'
    print(label_fit)
    yth = lorentzian(kx,popt[0],popt[1])
    
    fig, ax = plt.subplots()
    ax.plot(kx,y_exp,'o-')
    ax.plot(kx,yth,'r',label = label_fit)
    ax.legend()
    
    ax.set_xlabel(r'$k_x \; \mathrm{(rad.m^{-1})}$')
    ax.set_ylabel(r'$|\hat{V_x}|(k_x,k_y^{max}) \; \mathrm{(u.a.)}$')
    
    
    figname = f'{fig_folder}Lorentzian_fit_kx_{ID_txt}'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    
    
    plt.close('all')
    # Save f_max, kx0, alpha and y_exp in a dictionnary
    
    dict_res = {}
    dict_res['f_demod'] = f_max
    dict_res['f_ex'] = float(f_ex)
    dict_res['k'] = popt[0]
    dict_res['err_k'] = err_coeff[0]
    dict_res['alpha'] = popt[1]
    dict_res['err_alpha'] = err_coeff[1]
    dict_res['lorentz'] = {}
    dict_res['lorentz']['kx'] = kx
    dict_res['lorentz']['y_exp'] = y_exp
    dict_res['EXP'] = data['EXP']
    dict_res['ID'] = data['ID']
    
    file2save = f'{current_path}PIV_attenuation_results_{ID_txt}.pkl'
    with open(file2save, 'wb') as pf:
        pickle.dump(dict_res,pf)
        


##################################################
#%% ############# Study step by step #############
##################################################

h = 7.5 # frasil thickness 
date = '2024_07_11'
path2data = f'U:/Aurore_frasil/{date}_e_{h}mm_laser/matData/'

folderlist = glob.glob(f'{path2data}*')



















#%% Create a small movie 
time_interval = data['SCALE']['ft']
cmax = 4e-2
ani = make_movie(data['Vx'],data['t'],data['x'],data['y'],parula_map,cmax,time_interval = time_interval)

file2save = f'{fig_folder}Movie_Vx_{ID_txt}.mp4'
ani.save(file2save)
print('Animation saved')




#%% Make small movie of demodulated field 
N = 100 # nb frames for animation
phi = np.linspace(0,2*np.pi,N)
cmax = 1e-2

slow_down_ratio = 10
time_interval = round(slow_down_ratio/f_max/N * 1e3) # in milliseconds

ani = animation_complex_field(demod_field,phi,data['x'],data['y'],parula_map,cmax,time_interval = 1e3)

file2save = f'{fig_folder}Movie_demodulated_Vx_{ID_txt}.mp4'
ani.save(file2save)
print('Animation saved')










#%% Load laser profile for this frequency 

root = 'Y:/Banquise/Aurore/amplitude_vagues_laser/'

key_thickness = str(h)
thck_string = f'{h}mm/'
path2laser_data = f'{root}{thck_string}'

file2load = f'{path2laser_data}{ID}/{ID}.pkl'
print(file2load)

with open(file2load,'rb') as pf:
    Q = pickle.load(pf)

x = np.arange(0,Q.shape[0],step = 1)/data['SCALE']['facq_pix'] # distance x in m
t = np.arange(0,Q.shape[1],step = 1)*data['SCALE']['ft']

#%% Superpose Vx and Vz 
profile = np.mean(data['Vx'],axis = 1) # average Vx over y coordinate
profile_demod = np.mean(demod_field,axis = 1) # average demodulated field over y coordinate
demod_laser = np.mean(Q*np.exp(1j*2*np.pi*f_max*t),axis = -1)
phase_deriv = 1j*2*np.pi*f_max
Vz_laser = demod_laser*phase_deriv


frame = 0
fig, ax = plt.subplots()
ax.plot(x,Q[:,frame])
ax.plot(x,np.real(demod_laser))
# ax.tick_params(axis = 'y',labelcolor = 'tab:blue')

ax1 = ax.twinx()
# ax.plot(x,np.real(Vz_laser))
ax1.plot(data['x'],np.real(profile_demod),color = 'tab:green')








#%% Demodulate laser profile at excitation frequency 

demod_laser = np.mean(Q*np.exp(1j*2*np.pi*f_max*t),axis = -1)
N = 100
phase = np.linspace(0,2*np.pi,N)
laser_evol = np.zeros((len(demod_laser),N),dtype = complex)
for k in range(N):
    laser_evol[:,k] = demod_laser*np.exp(1j*phase[k])

# Show profile as a movie
fig, ax = plt.subplots(figsize = (12,9))
slow_down_ratio = 50
time_interval = round(slow_down_ratio/f_max/N * 1e3) # in milliseconds
nb_frames = len(phase)

ani = animation_profile(fig, ax, np.real(laser_evol), x, phase, nb_frames, time_interval)

# Save the created animation
file2save = f'animation_laser_demod_nb_frames{nb_frames}_sdr{slow_down_ratio}_{ID_txt}'
file2save = file2save.replace('.','p')
file2save = f'{fig_folder}{file2save}.mp4'
ani.save(file2save)
print('Animation saved')








#%% Superpose Vx and Vz 
profile = np.mean(data['Vx'],axis = 1) # average Vx over y coordinate
profile_demod = np.mean(demod_field,axis = 1) # average demodulated field over y coordinate
demod_laser = np.mean(Q*np.exp(1j*2*np.pi*f_max*t),axis = -1)
phase_deriv = 1j*2*np.pi*f_max
Vz_laser = demod_laser*phase_deriv

fig, ax = plt.subplots()
ax.plot(x,np.real(demod_laser))
# ax.plot(x,np.real(Vz_laser))
ax.plot(data['x'],np.real(profile_demod))

#%% Superpose raw Vx / demodulated Vx / laser profile



frame = 0
fig, ax = plt.subplots()
ax.plot(x,Q[:,frame])
ax.plot(data['x'],profile[:,frame])
ax.plot(data['x'],np.real(profile_demod))

















#%% Show demodulated laser and demodulated PIV

profile = np.mean(data['Vx'],axis = 1) # average Vx over y coordinate
profile_demod = np.mean(demod_field,axis = 1) # average demodulated field over y coordinate

PIV_evol = np.zeros((profile_demod.shape[0],N),dtype = complex)
for k in range(N):
    PIV_evol[:,k] = profile_demod*np.exp(1j*phase[k])
    
# Show profile as a movie
fig, ax = plt.subplots(figsize = (12,9))
slow_down_ratio = 50
time_interval = round(slow_down_ratio/f_max/N * 1e3) # in milliseconds
nb_frames = len(phase)


dict_demod = {}
dict_demod['laser'] = {'data':laser_evol,'x':x}
dict_demod['PIV'] = {'data':PIV_evol,'x':data['x']}

def animation_profile_laser_PIV(fig,ax,dict_demod,phi,nb_frames,time_interval):
    
    """ Create and return animation of a line plot using matplotlib.animation
    Inputs : - fig: matplotlib figure
             - ax: axis object
             - dict_demod : dictionnary containing several keys. Each key should contain following arguments :
                 - data: numpy array, data to plot #dim0 : space, #dim1 : time
                 - x: numpy array, x-axis (space)
             - phi: numpy array, (phase array)
             - nb_frames : int, number of frames to show
             - time_interval : time between two consecutive frames
             
    Output : - ani: matplotlib animation object"""
    
    
    j0 = 0 # initial time index 
    line = {}
    for key in dict_demod.keys():
        data = dict_demod[key]['data']
        x = dict_demod[key]['x']
        line[key] = ax.plot(x,data[:,j0],label = key)[0]
    
    ax.set_title(r'$phi =' + '{:.2f}'.format(phi[j0]) + r' \; \mathrm{(rad)}$')
    ax.set_xlabel(r'$x \; \mathrm{(mm)}$')
    ax.set_ylabel(r'$\xi \; \mathrm{(mm)}$')
    ax.legend()
    
    # small function to update the current figure 
    def update_profile_plot(frame):
        line.set_xdata(x)
        line.set_ydata(data[:,frame])
        ax.set_title(r'$phi =' + '{:.2f}'.format(phi[frame]) + r' \; \mathrm{(rad)}$')
        return line
    
    # create an animation 
    ani = animation.FuncAnimation(fig=fig, func=update_profile_plot, frames=nb_frames, interval=time_interval)
    plt.show()
    print('Animation computed')

    return ani







ani = animation_profile(fig, ax, np.real(laser_evol), x, phase, nb_frames, time_interval)

# Save the created animation
file2save = f'animation_laser_demod_nb_frames{nb_frames}_sdr{slow_down_ratio}_{ID_txt}'
file2save = file2save.replace('.','p')
file2save = f'{fig_folder}{file2save}.mp4'
ani.save(file2save)
print('Animation saved')
    









#%% Perform space time FFT

Efk = FT.space_time_spectrum(data['Vx'],1/data['SCALE']['fx'],data['SCALE']['facq_t'],add_pow2 = [1,1,0])

#%% Plot space-time spectrum 

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.imshow(Efk['E'], cmap = parula_map , aspect = 'auto', norm = 'log', 
              origin = 'lower', interpolation = 'gaussian',
              extent = (Efk['k'].min(),Efk['k'].max(),Efk['f'].min(),Efk['f'].max()))

ax.set_xlabel(r'$k \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$f \; \mathrm{(Hz)}$', labelpad = 5)

cbar = plt.colorbar(c,ax = ax)
cbar.set_label(r'$|\hat{V}_x| (k,\omega) \; \mathrm{(u.a.)}$',labelpad = 5)

# figname = f'{fig_folder}Efk_raw_{date}_{drone_ID}_{exp_ID}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.svg', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')




























