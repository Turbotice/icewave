# -*- coding: utf-8 -*-
"""
Created on Sun Jan 12 16:52:44 2025

@author: sebas
"""

import matplotlib.pyplot as plt
import matplotlib as mpl

#%% 

def set_matplotlib_param(style):
    
    # fig_size = (6.4,4.8) default fig_size 
    
    if style == 'single':
        fig_size = (8,6)
        font_size_medium = 22
    elif style == 'double':
        fig_size = (8,6)
        font_size_medium = 30
    elif style == 'triple':
        fig_size = (8,6)
        font_size_medium = 28
    elif style == 'square':
        fig_size = (12,9)
        font_size_medium = 20
    elif style == 'powerpoint':
        fig_size = (12,9)
        font_size_medium = 40
        mpl.rcParams['lines.markersize'] = 12
    elif type(style) == int or type(style) == float:
        font_size_medium = style
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
    plt.rc('font', family='serif', serif='Computer Modern')
    
