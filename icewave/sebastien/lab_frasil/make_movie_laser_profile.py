# -*- coding: utf-8 -*-
"""
Created on Fri May 23 14:59:46 2025

@author: sebas
"""

import matplotlib.pyplot as plt
import matplotlib.animation as animation

import pickle
import glob

from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# import icewave.sebastien.set_graphs as set_graphs

def set_matplotlib(style):
    # fig_size = (6.4,4.8) default fig_size 
    
    if style == 'single':
        fig_size = (8,6)
        font_size_medium = 22
    elif style == 'double':
        fig_size = (8,6)
        font_size_medium = 30
    elif style == 'square':
        fig_size = (12,9)
        font_size_medium = 20
    else:
        
        raise ValueError('The chosen style is not defined')
        
    plt.rcParams['figure.figsize'] = fig_size

    font_size_small = round(0.75*font_size_medium)
    plt.rc('font', size=font_size_medium)          # controls default text sizes
    plt.rc('axes', titlesize=font_size_medium)     # fontsize of the axes title
    plt.rc('axes', labelsize=font_size_medium)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=font_size_small)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=font_size_small)    # fontsize of the tick labels
    plt.rc('legend', fontsize=font_size_medium)    # legend fontsize
    plt.rc('figure', titlesize=font_size_medium)  # fontsize of the figure title

    plt.rcParams["axes.labelpad"] = 5
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='serif', serif='Computer Modern')
    return

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
    
    set_matplotlib('single')
    fig, ax = plt.subplots()
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

def process_file(args):
    
    file2load,path2data,nb_periods,slow_ratio = args
    
    with open(file2load,'rb') as pf:
        data = pickle.load(pf)

    h = data['h']
    f_ex = data['f_ex']
    amplitude = data['amplitude']

    suffixe = f'h_{h}_fex_{f_ex}_amp_{amplitude}mm'
    suffixe = suffixe.replace('.','p')

    nb_frames = int(nb_periods*data['SCALE']['facq_t']/data['f_ex'])
    ani = animation_profile(data['spatio']*1e3,data['x'],data['t'],nb_frames,
                            time_interval = 1e3*slow_ratio/data['SCALE']['facq_t'])

    file2save = f'{path2data}Laser_extraction/animation_{suffixe}.mp4'
    ani.save(file2save)
    
    return 

#%%
def main():
    h = 5.0 # frasil thickness 
    date = '2024_07_10'
    base = '/media/turbots/Backup25/'
    path2data = f'{base}Aurore_frasil/{date}_e_{h}mm_laser/'
    
    filelist = glob.glob(f'{path2data}Laser_extraction/scaled_laser_structure_*.pkl')
    
    nb_periods = 20
    slow_ratio = 5
    
    args_list = [(file,path2data,nb_periods,slow_ratio) for file in filelist]
    
    # for args in args_list:
    #     process_file(args)
        
    with ProcessPoolExecutor(max_workers=30) as executor:
        list(tqdm(executor.map(process_file, args_list), total=len(args_list)))

if __name__ == '__main__':
    main()

