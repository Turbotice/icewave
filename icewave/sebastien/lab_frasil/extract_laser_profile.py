# -*- coding: utf-8 -*-
"""
Created on Thu May  1 08:03:09 2025

@author: sebas
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation

import pickle
import os
import glob
import re
import cv2 as cv 

from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

# import icewave.tools.datafolders as df

# plt.rc('text', usetex=True)
# plt.rc('font', family='serif', serif='Computer Modern')

#%% FUNCTION SECTION 
def subpix_precision(profile, idx_max):
    """ Computes the subpixel precision of a position in the convolution product 
    cf thesis manuscript Antonin Marchand """
    
    p = profile[idx_max]
    p_right = profile[idx_max + 1]
    p_left = profile[idx_max - 1]

    delta = (p_right - p_left)/(2*(2*p - p_right - p_left))
    i_subpix = idx_max + delta # subpixel index 
    
    max_subpix = p + 0.5*(2*p - p_right - p_left)*delta**2 # max value at the subpixel idx_max
    
    return i_subpix,max_subpix

def img_correlation_along_axis(M,analysis,axis = 0):
    """ Compute correlation between each column of a 2D matrix and an analysis function. 
    An image is built from the analysis function (using np.tile) in order to perform correlation 
    
    Inputs : - M, numpy array, [ny,nx]
             - analysis, 1D numpy array, analysis function [N]
             - axis, int, to define on which axes correlations must be performed
             
    Outputs : - correlated, numpy array [ny,nx] containing convolution result"""
    
    reps = (M.shape[1],1)
    analysis_matrix = np.tile(analysis,reps)
    
    correlated = scipy.signal.fftconvolve(M,analysis_matrix.T,mode = 'same',axes = axis)
    return correlated
    
def extract_laser_frame(M,analysis):
    """ Extract laser for a single frame 
    Inputs : - M, numpy array [ny,nx], matrix from which laser position should be extracted
             - analysis, 1D numpy array, analysis function [N] 
             
    Outpouts : - laser_profile, numpy array shape (nx,), containing subpix positions of the laser """

    correlated = img_correlation_along_axis(M, analysis)
    
    laser_profile = np.zeros(correlated.shape[1])

    for j in range(correlated.shape[1]):
        # compute sub-pixel position 
        idx_max = np.argmax(correlated[:,j])
        idx_subpix,max_subpix = subpix_precision(correlated[:,j],idx_max)
        laser_profile[j] = idx_subpix
        
    return laser_profile

def get_laser_indices(img,threshold,channel_color):
    """ Get idx_start and idx_end, indices of columns over which we are looking for the laser
    Inputs : - img, image [ny,nx,3]
             - threshold, float, threshold used to binarize image
             - channel_color, 0,1 or 2, channel over which we perform thresholding 0 = R, 1 = G, 2 = B 
    Outputs : idx_start and idx_end, indices of column between which we are looking for the laser """
    
    mask = img[:,:,channel_color] > threshold

    # compute sum along raws
    sum_mask = np.sum(mask,axis = 0)
    # get minimal and maximal index between which we perform correlation 
    laser_indices = np.where(sum_mask > 0)[0]
    idx_start = laser_indices[0]
    idx_end = laser_indices[-1]
    return idx_start,idx_end


def scale_spatio(laser_spatio,theta_rad,facq_pix,facq_t,idx_start,idx_end):
    """ Scale spatio temporal extraction of laser 
    Inputs : - laser_spatio, numpy array [nx,nt], spatio temporal extraction of the laser
             - theta, float, laser angle to the vertical, in radians
             - facq_pix, float, spatial acquisition frequency, in pix/meter
             - facq_t, float, temporal acquisition frequency, in frame/seconds
             - idx_start and idx_end, column indices of the image between which we looked for the laser 
             
    Outputs : - data, dictionnary containing parameters and scalings but mainly the 3 main keys :
                  + spatio_scaled, numpy array [nx,nt], surface deformation spatio temporal
                  + x, 1D numpy array, x-coordinate
                  + t, 1D numpy array, time coordinate """

    # scaling 
    spatio_scaled = (laser_spatio - np.mean(laser_spatio))/np.tan(theta_rad)/facq_pix # deformation in meter 
    x = np.arange(spatio_scaled.shape[0])/facq_pix # x coordinate in meter
    t = np.arange(spatio_scaled.shape[1])/facq_t # time array in seconds 
    
    # structure data 
    data = {}
    data['spatio'] = spatio_scaled
    data['x'] = x
    data['t'] = t
    data['param'] = {}
    data['param']['idx_start'] = idx_start
    data['param']['idx_end'] = idx_end
    data['param']['theta'] = theta_rad
    
    data['SCALE'] = {}
    data['SCALE']['facq_t'] = facq_t
    data['SCALE']['facq_pix'] = facq_pix
    data['SCALE']['fx'] = 1/facq_pix
    data['SCALE']['units'] = {'facq_t':'frame/s','facq_pix':'pix/m','fx':'m/pix'}
    
    return data 

def single_experiment_extraction(path2img,table_param,results_folder):
    """ perform laser extraction for all frames of a single movie 
    Inputs :- path2img, string, path to all frames 
            - table_param, dictionnary containing all experimental parameters
            - results_folder, string, folder where results will be saved """
            
    # generate a suffixe for saving data 
    h = float(re.findall(r'e_(\d+\.\d+)mm',path2img)[0]) # frasil thickness
    f_ex = float(re.findall(r'(\d+\.\d+)Hz',path2img)[0]) # excitation_frequency 
    amplitude = float(re.findall(r'Hz_(\d+)mm',path2img)[0])
    
    # collect parameters for this experiment
    key = f'h_{h}mm_fex_{f_ex}Hz'
    theta = table_param[key]['theta']
    theta_rad = theta*np.pi/180
    facq_t = table_param[key]['facq_t']
    facq_pix = table_param[key]['facq_pix']
    
    
    suffixe = f'h_{h}_fex_{f_ex}_amp_{amplitude}mm'
    suffixe = suffixe.replace('.','p')
    
    # list of image 
    img_list = glob.glob(f'{path2img}Basler*.tiff')

    ###### Define section of image where laser will be extracted
    red_threshold = 120 # threshold for selection of idx_start and idx_end
    channel_color = 0
    file2load = img_list[0]
    img = cv.imread(file2load)
    img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
    idx_start,idx_end = get_laser_indices(img, red_threshold, channel_color)
    
    fig, ax = plt.subplots()
    ax.imshow(img)
    ax.axvline(idx_start)
    ax.axvline(idx_end)
    figname = f'{results_folder}Subimage_selection_{suffixe}'
    plt.savefig(figname + '.pdf', bbox_inches='tight')
    plt.savefig(figname + '.svg', bbox_inches='tight')
    plt.savefig(figname + '.png', bbox_inches='tight')
    
    plt.close(fig)
    
    ####### Loop over all images 
    laser_spatio = np.zeros((idx_end - idx_start,len(img_list)))
    
    gaussian_width = 15
    N = 10
    gaussian = scipy.signal.windows.gaussian(M = gaussian_width * N,std = gaussian_width)
    
    for idx_img in range(len(img_list)):
        
        file2load = img_list[idx_img]

        img = cv.imread(file2load)
        img = cv.cvtColor(img,cv.COLOR_BGR2RGB)
            
        img_laser = img[:,idx_start:idx_end,:]
    
        laser_profile = extract_laser_frame(img_laser[:,:,0],gaussian)
        laser_spatio[:,idx_img] = laser_profile

    file2save = f'{results_folder}unscaled_laser_{suffixe}.pkl'

    with open(file2save,'wb') as pf:
        pickle.dump(laser_spatio,pf)
        

    file2save = f'{results_folder}unscaled_laser_{suffixe}.pkl'
    with open(file2save,'rb') as pf:
        laser_spatio = pickle.load(pf)
        
    data = scale_spatio(laser_spatio,theta_rad,facq_pix,facq_t,idx_start,idx_end)

    data['h'] = h
    data['f_ex'] = f_ex    
    data['amplitude'] = amplitude

    file2save = f'{results_folder}scaled_laser_structure_{suffixe}.pkl'
    with open(file2save,'wb') as pf:
        pickle.dump(data,pf)
        
    return 

def process_folder(args):
    folder, table_param, results_folder = args
    path2img = f'{folder}/'
    single_experiment_extraction(path2img, table_param, results_folder)
    print(f'Laser extraction performed for {path2img}')

#%% Load image
def main():
    
    h = 10.0 # frasil thickness 
    date = '2024_07_11'
    # base = df.find_path(disk='Backup25',year='2024')
    base = '/media/turbots/Backup25/'
    path2data = f'{base}Aurore_frasil/{date}_e_{h}mm_laser/'
    
    folderlist = glob.glob(f'{path2data}*mm')
    
    print(folderlist)
    # create results folder 
    results_folder = f'{path2data}Laser_extraction/'
    if not os.path.isdir(results_folder) :
        os.mkdir(results_folder)
        
    # load experimental parameters dictionnary
    path2table = f'{base}Aurore_frasil/table_experimental_parameters.pkl'
    with open(path2table,'rb') as pf :
        table_param = pickle.load(pf)
           
    # select an experiment (an excitation frequency)
    # for folder in tqdm(folderlist):
    #     process_folder(folder)
    
    # create a list of arguments
    args_list = [(folder, table_param, results_folder) for folder in folderlist]

    with ProcessPoolExecutor(max_workers=30) as executor:
        list(tqdm(executor.map(process_folder, args_list), total=len(args_list)))

if __name__ == '__main__':
    main()





