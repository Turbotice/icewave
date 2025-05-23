# -*- coding: utf-8 -*-
"""
Created on Fri May 23 14:59:46 2025

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import pickle
import os
import glob

import icewave.tools.matlab2python as mat2py
import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs
import icewave.tools.Fourier_tools as FT


def animation_profile(data,x,t,nb_frames,time_interval):
    
    """ Create and return animation of a line plot using matplotlib.animation
    Inputs : - fig: matplotlib figure
         - ax: axis object
         - data: numpy array, data to plot #dim0 : space, #dim1 : time
         - x: numpy array, x-axis (space)
         - t: numpy array, (time)
         - nb_frames : int, number of frames to show
         - time_interval : time between two consecutive frames
         
    Output : - ani: matplotlib animation object"""
    
    set_graphs.set_matplotlib_param('single')
    fig, ax = plt.subplots(figsize = (12,9))
    j0 = 0 # initial time index 
    line = ax.plot(x,data[:,j0])[0]
    
    ax.set_title(r'$t =' + '{:.2f}'.format(t[j0]) + r' \; \mathrm{s}$')
    ax.set_xlabel(r'$x \; \mathrm{(m)}$',labelpad = 5)
    ax.set_ylabel(r'$\xi \; \mathrm{(mm)}$',labelpad = 5)
    
    # small function to update the current figure 
    def update_profile_plot(frame):
        line.set_xdata(x)
        line.set_ydata(data[:,frame])
        ax.set_title(r'$t =' + '{:.2f}'.format(t[frame]) + r' \; \mathrm{s}$')
        return line
    
    # create an animation 
    ani = animation.FuncAnimation(fig=fig, func=update_profile_plot, frames=nb_frames, interval=time_interval)
    plt.show()
    print('Animation computed')
    
    return ani

#%%

h = 10.0 # frasil thickness 
date = '2024_07_11'
path2data = f'U:/Aurore_frasil/{date}_e_{h}mm_laser/'

filelist = glob.glob(f'{path2data}Laser_extraction/scaled_laser_structure_*.pkl')

idx_experiment = 5
file2load = filelist[idx_experiment]
print(file2load)
with open(file2load,'rb') as pf:
    data = pickle.load(pf)

#%% Compute animation and save it 

nb_periods = 20
slow_ratio = 5
nb_frames = int(nb_periods*data['SCALE']['facq_t']/data['f_ex'])
ani = animation_profile(data['spatio']*1e3,data['x'],data['t'],nb_frames,
                        time_interval = 1e3*slow_ratio/data['SCALE']['facq_t'])

file2save = f'{path2data}Laser_extraction/animation_test.mp4'
print('Saving animation...')
ani.save(file2save)
print('Animation saved')