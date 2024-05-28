# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 11:14:08 2023

@author: sebas
"""


import matplotlib.pyplot as plt

import glob
import numpy as np
import cv2 as cv
import os
import pickle
import csv
### Enable Latex in Matplotlib
plt.rcParams.update({
    "text.usetex": True,
})

from scipy.optimize import curve_fit
from scipy import signal
from scipy.fft import fft2, fftshift, fftfreq
from scipy.signal import argrelextrema
from skimage.measure import label, regionprops
from scipy.optimize import curve_fit

import seb
from seb.pickle_m import read, write

#%% Basics functions 

def detect_fronts(img_particle, analyzing_func, col, x_crop) :
    """ Detect and return the fronts of the bottom of the particle along a given column
    
    It takes as arguments : 
        - img_particle, the image of the float in gray scale 
        - analyzing_func, a 1D array that represents the discretized analyzing function
        - col, the column along which the fronts are detected 
        - x_crop, the boundaries of the column 
        
    It returns a tuple of the detected back and front edges of the object """
    
    (x_crop_min, x_crop_max) = x_crop
    
    # crop the image and get the grey profile along the column col
    grey_profile = img_particle[x_crop_min:x_crop_max, col]
    # convolve an analyzing function to the grey profile
    convolution = np.convolve(grey_profile,analyzing_func, 'same')
    
    # detect the position of the two maximas of the convolution product
    local_max = argrelextrema(convolution,np.greater)
    max_corr_1 = int(np.where(convolution == np.sort(convolution[local_max])[-1])[0])
    max_corr_2 = int(np.where(convolution == np.sort(convolution[local_max])[-2])[0])

    max_corr_1 = subpix_precision(convolution, max_corr_1)
    max_corr_2 = subpix_precision(convolution, max_corr_2)
    x_back = min(max_corr_1, max_corr_2) + x_crop_min
    x_front = max(max_corr_1, max_corr_2) + x_crop_min
    
    return (x_back,x_front)
    
def subpix_precision(convolution, idx_max):
    """ Computes the subpixel precision of a position in the convolution product 
    cf thèse de Antonin Marchand """
    
    p = convolution[idx_max]
    p_right = convolution[idx_max + 1]
    p_left = convolution[idx_max - 1]

    delta = (p_right - p_left)/(2*(2*p - p_right - p_left))
    i_subpix = idx_max + delta
    
    return i_subpix

def keep_coherent_pts(front_array, column_array, min_difference = 10):
    """ Function that search for incoherent points and get rid off them
    A point is said to be incoherent when distance (in pixels) from the polynomial fit is larger 
    than the min_difference"""

    idx = np.arange(0,len(front_array))
    # approximation by a degree 2 polynomial
    coef = np.polyfit(idx,front_array,2)
    poly = coef[2] + coef[1]*idx + coef[0]*idx**2

    # detect points that are too far from the fit :
    diff = np.abs(front_array - poly)
    outliers = np.where(diff > min_difference)[0]

    # get rid of outliers
    new_front = np.delete(front_array, outliers)
    new_columns = np.delete(column_array, outliers)
    return new_front, new_columns

def y_displacement(a, alpha = -3.159, beta = 882.75):
    """ Compute the displacement along y(t), in mm, from the evolution of the semi-major axis
    
    Alpha and beta were calculated from the code 'Calibration' """
    y = alpha*a + beta
    return y

def reducing_signal(arr):
    """ Reduce the mean value from the entered signal and the linear part of the signal """
    s = arr - np.mean(arr)
    idx = np.arange(0,len(s))
    # approximation by a 1 degree polynomial
    coef = np.polyfit(idx,s,1)
    poly = coef[1] + coef[0]*idx
    return s - poly

def write_frequencies(save_path, suffixe, signal_label, main_peaks_propreties):
    """ Write main frequencies from a Fourier spectrum in a txt file 
    - save_path : path at which the file is saved 
    - suffixe : used to name the txt file
    - signal_label : name of the variable we are looking at 
    - main_peaks propreties : (main_freq, main_amplitudes) of the main peaks of the spectrum """
    
    # create text file
    text_file = save_path + 'frequencies_' + suffixe +'.txt'
    main_freq, main_amplitude = main_peaks_propreties
    # write heading 
    with open(text_file, 'w') as f:
        f.write('Main frequencies of ' + signal_label)
        f.write('\n')
        f.write('Frequency'+ ' ' +'Amplitude')
        f.write('\n')
    f.close()   
    # write frequency and amplitude of the peaks
    with open(text_file, 'a') as f:
        for i in range(len(main_freq)-1,-1,-1):
            f.write(str("{:.3f}".format(main_freq[i])) + '  ' + str("{:.5f}".format(main_amplitude[i])))
            f.write('\n')
    f.close()
    
def write_frequencies_csv(save_path, suffixe, signal_label, main_peaks_propreties):
    """ Write main frequencies from a Fourier spectrum in a txt file 
    - save_path : path at which the file is saved 
    - suffixe : used to name the txt file
    - signal_label : name of the variable we are looking at 
    - main_peaks propreties : (main_freq, main_amplitudes) of the main peaks of the spectrum """
    
    # create text file
    text_file = save_path + 'frequencies_' + suffixe +'.csv'
    main_freq, main_amplitude = main_peaks_propreties
    # write heading 
    with open(text_file, 'w', newline = '') as f:
        writer = csv.writer(f)
        writer.writerow(['Frequency', 'Amplitude'])  
    # write frequency and amplitude of the peaks
        for i in range(len(main_freq)-1,-1,-1):
            writer.writerow([main_freq[i], main_amplitude[i]])
    f.close()

def Fourier_spectrum(s, f_sampling, signal_label, save_path = '', suffixe = '', h = 0):
    """ Computes and plot the Fourier spectrum of a signal s, using the following arguments :
        - signal s
        - f_sampling = fps/step, the frequency at which the signal is sampled 
        - h is the height of the main peaks that will be tagged (should be a tuple)
        - signal_label, the label of the features we are Fourier transforming """
        
    # zero padding 
    N = len(s)
    power_2 = int(np.log(N)/np.log(2))
    padding_length = 2**(power_2 + 3)

    # computes fourier transform with zero padding
    fourier = np.fft.fft(s, n = padding_length, norm = 'forward')
    fourier_abs = np.abs(fourier[:padding_length//2])
    # computes associated frequencies
    freq = np.fft.fftfreq(padding_length, d = 1/f_sampling)
    # keep only positive frequencies 
    freq_pos = freq[:padding_length//2]
    #keep only frequencies below 10 Hz
    f_plot = freq_pos[freq_pos < 10]
    fourier_plot = fourier_abs[:len(f_plot)]
    
    fig, ax = plt.subplots()
    ax.plot(f_plot,fourier_plot)
    
    ### Find peaks of the Fourier transform
    # find peaks higher than h, separated by more than distance samples from each other
    peaks_index, peaks_propreties = signal.find_peaks(fourier_plot, height = h, distance = 1)
    peaks_amplitude = peaks_propreties['peak_heights']
    peaks_freq = f_plot[peaks_index]
    # get the peaks with the bigger amplitude
    main_idx = np.argsort(peaks_amplitude)
    main_freq = peaks_freq[main_idx]
    main_amplitude = peaks_amplitude[main_idx]
    # plot main peaks 
    ax.plot(main_freq,main_amplitude,'ok')
    
    # annotate the graph
#     for i in range(len(main_idx)):
#         xy = (main_freq[i],main_amplitude[i])
#         xytext = (main_freq[i] + 0.5,main_amplitude[i])
#         ax.annotate(''+str("{:.3f}".format(main_freq[i]))+'',xy, xytext)

    ax.set_ylabel(r'Amplitude')
    ax.set_xlabel(r'Frequency (Hz)')
    ax.set_title('FFT of ' + str(signal_label) +'')
    ax.grid()
    
    # keep main peaks propreties
    main_peaks_propreties = (main_freq, main_amplitude)
    # write peaks propreties
    if save_path != '':
        write_frequencies(save_path, suffixe, signal_label, main_peaks_propreties)
        write_frequencies_csv(save_path, suffixe, signal_label, main_peaks_propreties)
    return ax

#%%   Functions used to fit an ellipse to a serie of points


def fit_ellipse(x, y):
    """

    Fit the coefficients a,b,c,d,e,f, representing an ellipse described by
    the formula F(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0 to the provided
    arrays of data points x=[x1, x2, ..., xn] and y=[y1, y2, ..., yn].

    Based on the algorithm of Halir and Flusser, "Numerically stable direct
    least squares fitting of ellipses'.


    """

    D1 = np.vstack([x**2, x*y, y**2]).T
    D2 = np.vstack([x, y, np.ones(len(x))]).T
    S1 = D1.T @ D1
    S2 = D1.T @ D2
    S3 = D2.T @ D2
    T = -np.linalg.inv(S3) @ S2.T
    M = S1 + S2 @ T
    C = np.array(((0, 0, 2), (0, -1, 0), (2, 0, 0)), dtype=float)
    M = np.linalg.inv(C) @ M
    eigval, eigvec = np.linalg.eig(M)
    con = 4 * eigvec[0]* eigvec[2] - eigvec[1]**2
    ak = eigvec[:, np.nonzero(con > 0)[0]]
    return np.concatenate((ak, T @ ak)).ravel()


def cart_to_pol(coeffs):
    """

    Convert the cartesian conic coefficients, (a, b, c, d, e, f), to the
    ellipse parameters, where F(x, y) = ax^2 + bxy + cy^2 + dx + ey + f = 0.
    The returned parameters are x0, y0, ap, bp, e, phi, where (x0, y0) is the
    ellipse centre; (ap, bp) are the semi-major and semi-minor axes,
    respectively; e is the eccentricity; and phi is the rotation of the semi-
    major axis from the x-axis in rad

    """

    # We use the formulas from https://mathworld.wolfram.com/Ellipse.html
    # which assumes a cartesian form ax^2 + 2bxy + cy^2 + 2dx + 2fy + g = 0.
    # Therefore, rename and scale b, d and f appropriately.
    a = coeffs[0]
    b = coeffs[1] / 2
    c = coeffs[2]
    d = coeffs[3] / 2
    f = coeffs[4] / 2
    g = coeffs[5]

    den = b**2 - a*c
    if den > 0:
        raise ValueError('coeffs do not represent an ellipse: b^2 - 4ac must'
                         ' be negative!')

    # The location of the ellipse centre.
    x0, y0 = (c*d - b*f) / den, (a*f - b*d) / den

    num = 2 * (a*f**2 + c*d**2 + g*b**2 - 2*b*d*f - a*c*g)
    fac = np.sqrt((a - c)**2 + 4*b**2)
    # The semi-major and semi-minor axis lengths (these are not sorted).
    ap = np.sqrt(num / den / (fac - a - c))
    bp = np.sqrt(num / den / (-fac - a - c))

    # Sort the semi-major and semi-minor axis lengths but keep track of
    # the original relative magnitudes of width and height.
    width_gt_height = True
    if ap < bp:
        width_gt_height = False
        ap, bp = bp, ap

    # The eccentricity.
    r = (bp/ap)**2
    if r > 1:
        r = 1/r
    e = np.sqrt(1 - r)

    # The angle of anticlockwise rotation of the major-axis from th x-axis (vertical here), in rad
    if b == 0:
        phi = 0 if a < c else np.pi/2
    else:
        phi = np.arctan((2.*b) / (a - c)) / 2
        if a > c:
            phi += np.pi/2
    if not width_gt_height:
        # Ensure that phi is the angle to rotate to the semi-major axis.
        phi += np.pi/2
    phi = phi % np.pi


    return x0, y0, ap, bp, e, phi


def get_ellipse_pts(params, npts=100, tmin=0, tmax=2*np.pi):
    """
    Return npts points on the ellipse described by the params = x0, y0, ap,
    bp, e, phi for values of the parametric variable t between tmin and tmax.

    """

    x0, y0, ap, bp, e, phi = params
    # A grid of the parametric variable, t.
    t = np.linspace(tmin, tmax, npts)
    x = x0 + ap * np.cos(t) * np.cos(phi) - bp * np.sin(t) * np.sin(phi)
    y = y0 + ap * np.cos(t) * np.sin(phi) + bp * np.sin(t) * np.cos(phi)
    return x, y


#%%  MAIN FUNCTIONS TO PROCESS A VIDEO

def process_frame(frame, analyzing_func, min_area = 15000 , boarder = 40, peaks_separation = 20):
    """ This function find the float on the frame and fit the bottom of the float with an ellipse, it takes as arguments :
    - frame, an image in grayscale
    - analyzing_func, the analyzing function used for the fronts detection
    - min_area, minimal pixel area to detect the particle 
    - boarder, default, the boarder added to the image of the particle
    - peaks_separation, default, the distance (in pixels) from the float sides at which the fronts detection along columns
    is performed 
    
    It returns : coordinates of the ellipse center (in the camera plane),
    its major and minor semi axis and its angle with the vertical (x axis) """
    
    ### Identification of the float
    ret,img_otsu = cv.threshold(frame,0, 255, cv.THRESH_BINARY_INV+cv.THRESH_OTSU) # optimized thresholding using Otsu's method

    label_image = label(img_otsu)
    for region in regionprops(label_image):
        if region.area>min_area :
            circle = region
            print(region.area)

    xmin = np.min(circle.coords[:,0])
    xmax = np.max(circle.coords[:,0])
    ymin = np.min(circle.coords[:,1])
    ymax = np.max(circle.coords[:,1])
    
    # crop the frame using the tilt of the object
 
    # get the column idx associated to xmin
    col_xmin = np.mean(circle.coords[np.where(circle.coords[:,0] == xmin)[0],1])
    # get the closest distance to ymin and ymax
    distance = np.asarray([abs(col_xmin - ymin), abs(col_xmin - ymax)])
    closest = np.argmin(distance)
    # closest = 0 if col_xmin is closer to ymin, 1 if it is closer to ymax
    
    if distance[closest] < 50 :
        if closest : # closest to ymax
            # get row associated to ymax (right boarder of props)
            row_ymax = circle.coords[np.where(circle.coords[:,1] == ymax)[0],0][-1]
            # if the object is too tilted, resize the img_particle
            if abs(row_ymax - xmin) > 20 :
                ymax = circle.coords[np.where(circle.coords[:,0] == xmin)[0],1][-1] + 20
        else : # closest to ymin
            # get row associated to ymin (left boarder of props)
            row_ymin = circle.coords[np.where(circle.coords[:,1] == ymin)[0],0][0]
            # if the object is too tilted, resize the img_particle
            if abs(row_ymin - xmin) > 20 :
                ymin = circle.coords[np.where(circle.coords[:,0] == xmin)[0],1][0]


    img_particle = frame[xmin - boarder:xmax + boarder,ymin - boarder:ymax + boarder]
    
    ### Fronts detection 
    # array of columns through which peaks detection is done
    column_array = np.arange(boarder + peaks_separation, ymax - ymin + boarder - peaks_separation)
    # array giving the position of back and fronts positions
    # row 0 corresponds to back 
    # row 1 corresponds to front 
    fronts_array = np.zeros((2,len(column_array)))
    
    x_crop = (0, boarder + (xmax - xmin)//2) # define the rows on which the front detection is performed
    for i in range(len(column_array)) :
        (fronts_array[0,i],fronts_array[1,i]) = detect_fronts(img_particle, analyzing_func, column_array[i], x_crop)
    
    ### Fit the bottom of the float with an ellipse
    # keep only coherent points
    back, column_back = keep_coherent_pts(fronts_array[0,:], column_array)
    front, column_front = keep_coherent_pts(fronts_array[1,:], column_array)

    ### get x and y coordinates used to build the ellipse
    x_ellipse = np.concatenate((back,front))
    y_ellipse = np.concatenate((column_back,column_front))

    # get coefficients of the ellipse :
    coeffs = fit_ellipse(x_ellipse,y_ellipse)
    x0, y0, ap, bp, e, phi = cart_to_pol(coeffs)
    
    # get position of the ellipse in the real frame
    x_center = xmin - boarder + x0
    y_center = ymin - boarder + y0
    
#     # get angle to the horizontal (clockwise)
#     theta = np.pi - phi
    return x_center, y_center, ap, bp, phi 


def process_tracking(video, save_path, analyzing_func, fps, variables, min_area, boundary_frames = (0,-1), x_crop = (0,-1), y_crop = (0,-1)):
    """ Process a sequence of images, for each frame, the bottom of the floatter is tracked and fit by an ellipse 
    Arguments :
    - paths[0] : video_path, paths[1] : save_path
    - analyzing_func : function used to detect fronts
    - fps : frame rate used
    - variables : dictionnary of variables 
    
    The function returns the features of the float bottom with time : 
    - x_center and y_center, the coordinates of the ellipse center 
    - ap and bp, the major and minor semi-axis of the ellipse
    - phi, the angle (in rad) of the ellipse with the horizontal """
    
#     video_path = path[0]
#     save_path = path[1]
    
#     if os.path.isdir(save_path) == False:
#         os.mkdir(save_path)
    
    ### Loads all the image 
#     video = glob.glob(video_path)
    video.sort()
    # create a features array to follow characteristics of the float with time
    first_frame, last_frame = boundary_frames
    video = video[first_frame : last_frame]
    features = np.zeros((5,len(video)))
    
    (x_crop_min, x_crop_max) = x_crop
    (y_crop_min, y_crop_max) = y_crop
    for i in range(len(video)) :
        frame = cv.imread(video[i],cv.IMREAD_GRAYSCALE)
        frame = frame[x_crop_min:x_crop_max,y_crop_min:y_crop_max]
        print(i)
        features[:,i] = process_frame(frame, analyzing_func, min_area)
    
    ## smooth the semi-major axis 
    diameter = variables['d']
    # apply a Savitzky - Golay filter to the evolution of the semi-major axis 
    smooth_a = signal.savgol_filter(features[2,:],51,3)
    # create the scale factor in mm/pixels for each frame
    fx = diameter/smooth_a/2
    # computes y displacement
    y = y_displacement(features[2,:])
    #
    features_dict = {'x':features[1,:]*fx, 'y' : y, 'z':features[0,:]*fx,
                    'a': features[2,:]*fx, 'b': features[3,:]*fx, 'theta': np.pi/2 - features[4,:]}
    # create suffixe
    chain = 'data_'
    for key in variables.keys():
        chain += key + str(round(variables[key],3)) + '_'
        suffixe = chain.replace('.','p')
        
    # save the dictionnary as a pickle 
    save_file = save_path + suffixe + '.pkl'
    write(features_dict,save_file)
    return features_dict

def plot_features(features_dict, save_path, fps, variables, scale = 0):
    """ Function that plots and saves features evolution of the ellipse that matches the float bottom 
    - features_dict : different values that we track and evaluate (dictionnary)
    - fps : frame rate used on the video 
    - save_path : path at which we want to save our plots
    - variables : a dictionnary with the parameters of the experiment"""
    
    exp_caract = 'param :'
    for key in variables.keys():
        exp_caract += key + str(round(variables[key],3)) + ', '
        
    # create an array of time
    time = np.arange(0,len(features_dict['x']))/fps
    
    # plot evolution of the position 
    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(time,features_dict['x'],'-',label = 'x(t)')
    ax.set_title('Displacement x(t), '+ exp_caract)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('x (mm)')
    ax.grid()
    
    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(time,features_dict['y'],'-',label = 'y(t)')
    ax.set_title('Displacement y(t), '+ exp_caract)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('y (mm)')
    ax.grid()   

    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(time,features_dict['z'],'-',label = 'z(t)')
    ax.set_title('Displacement z(t), '+ exp_caract)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('z (mm)')
    ax.grid()
    
    # plot evolution of the semi-axis
    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(time,features_dict['a'],'-',label = 'a(t)')
    ax.set_title('Semi-major axis a(t), '+ exp_caract)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('a (mm)')
    ax.grid()

    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(time,features_dict['b'],'-',label = 'b(t)')
    ax.set_title('Semi-minor axis b(t), '+ exp_caract)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('b (mm)')
    ax.grid()
    
    # plot evolution of the ellipse orientation
    fig, ax = plt.subplots(figsize = (10,10))
    ax.plot(time,features_dict['theta'],'-',label = 'theta(t)')
    ax.set_title(r'Angle $\theta$(t), '+ exp_caract)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel(r'$\theta$ (rad)')
    ax.grid()

def Fourier_features(features_dict, f_sampling, padding_length = 2**15):
    """ Return the FFT of each features of the features_dict. The datas are sampled with a frequency f_sampling
    - features_dict : dictionnary of features
    - f_sampling : frequency at which features are recorded """
    
    fft_dict = {}
    
    for key in features_dict.keys():
        s = reducing_signal(features_dict[key])
        N = len(s)
        # zero padding 
        
        if padding_length < N :
            power_2 = int(np.log(N)/np.log(2))
            padding_length = 2**(power_2 + 3)
            print('Padding length is too small for this signal')

        # computes fourier transform with zero padding
        fourier = np.fft.fft(s, n = padding_length, norm = 'forward')
        fourier_abs = np.abs(fourier[:padding_length//2])
        # computes associated frequencies
        freq = np.fft.fftfreq(padding_length, d = 1/f_sampling)
        # keep only positive frequencies 
        freq_pos = freq[:padding_length//2]
        #keep only frequencies below 10 Hz
        f_plot = freq_pos[freq_pos < 10]
        fourier_plot = fourier_abs[:len(f_plot)]
        
        # add fourier value to the dictionnary
        fft_dict[key] = fourier_plot
    
    # add the array of frequencies to the dictionnary
    fft_dict['freq'] = f_plot
    
    return fft_dict

def FFT_z(features_dict, fps, save_path, variables):
    """ Plot the FFT of the z displacement and keep frequencies in a .txt file 
    It takes as arguments : 
    - a features dictionnary, containing the z displacement 
    - fps, the frame rate used in our experiments 
    - save_path, the path where we can save the plot and the main frequencies 
    - variables, parameters used in our system """
    
    # get z displacement
    z = features_dict['z']
    signal = reducing_signal(z)
    
    # create filename to save
    chain = 'FFT_z_'
    for key in variables.keys():
        chain += key + str(round(variables[key],3)) + '_'
        suffixe = chain.replace('.','p')
    save_file = save_path + suffixe + '.pkl'
    
    # computes and plot the Fourier spectrum
    signal_label = suffixe.strip('FFT_z_')
    ax = Fourier_spectrum(signal, fps, signal_label, save_path, suffixe)
    # save the graph
    write(ax,save_file)
    pdf_file = save_file.strip('pkl')  + 'pdf'
    plt.savefig(pdf_file)
    
def variables_to_chain(variables_dict, prefixe):
    
    """" Transforms a dictionnary of parameters used for one experiment into a chain of character """
    
    chain = prefixe 
    for key in variables_dict.keys():
        chain += '_' + key + str(round(variables_dict[key],3)) 
    return chain
    
#%% Plot section 
### PLOTS FOR THE INTERNSHIP REPORT ### 

def plot_signal_and_Fourier(s,f_sampling,fig_path, title, signal_label, height = 0.005, padding_length = 2**15, text_size = 18):
    """ Plot a signal s and its Fourier transform in a single subplot and saves it
    
    It takes as arguments :
        - s the signal to plot
        - f_sampling, the sampling with wich the signal is sampled
        - fig_path, the path where the picture will be saved at format .jpg, .pdf 
        - title, the title of the graph (in latex)
        - height, peaks amplitude above which the peak of the Fourier transform will be labeled """


    plt.rcParams.update({'font.size': text_size})
    time = np.arange(0,len(s))/f_sampling
     
    s_copy = np.copy(s)
    s_copy = reducing_signal(s_copy)
    N = len(s)
    
    # zero padding 
    if padding_length < N :
        power_2 = int(np.log(N)/np.log(2))
        padding_length = 2**(power_2 + 3)
        print('Padding length is too small for this signal')
        
    # computes fourier transform with zero padding
    fourier = np.fft.fft(s_copy, n = padding_length, norm = 'forward')
    fourier_abs = np.abs(fourier[:padding_length//2])
    # computes associated frequencies
    freq = np.fft.fftfreq(padding_length, d = 1/f_sampling)
    # keep only positive frequencies 
    freq_pos = freq[:padding_length//2]
    #keep only frequencies below 10 Hz
    f_plot = freq_pos[freq_pos < 10]
    fourier_plot = fourier_abs[:len(f_plot)]
        
    ### Find peaks of the Fourier transform
    # find peaks higher than h, separated by more than distance samples from each other
    peaks_index, peaks_propreties = signal.find_peaks(fourier_plot, height, distance = 1)
    # peaks_amplitude = peaks_propreties['peak_heights']
    # peaks_freq = f_plot[peaks_index]
    # # get the peaks with the bigger amplitude
    # main_idx = np.argsort(peaks_amplitude)
    # main_freq = peaks_freq[main_idx]
    # main_amplitude = peaks_amplitude[main_idx]
    
    main_freq = np.zeros(len(peaks_index))
    main_amplitude = np.zeros(len(peaks_index))
    
    fig, ax = plt.subplots(1,2, figsize = (14,8))
    ax[0].plot(time,s)
    ax[0].set_xlabel(r'Time (s)')
    ax[0].set_ylabel(signal_label)
    ax[0].set_title(title)
    ax[0].tick_params(axis='both', which='major', labelsize = text_size - 2)
    
    ax[1].plot(f_plot,fourier_plot)
    ax[1].set_xlabel(r'Frequency (Hz)')
    ax[1].set_ylabel(r'Amplitude')
    ax[1].set_title(r'(b) FFT of the signal')
    ax[1].tick_params(axis='both', which='major', labelsize = text_size - 2)
    
    #annotate the graph
    for i in range(len(peaks_index)):
        # subpix precision
        idx = peaks_index[i]
        p = fourier_plot[idx]
        p_right = fourier_plot[idx+1]
        p_left = fourier_plot[idx-1]
    
        delta = (p_right - p_left)/(2*(2*p - p_right - p_left))
        new_freq = f_plot[idx] + (f_plot[idx+1] - f_plot[idx])*delta
        
        amplitude_subpix = p + 0.5*(2*p - p_right - p_left)*delta**2
        # keep the main frequencies and amplitude
        main_freq[i] = new_freq
        main_amplitude[i] = amplitude_subpix
        # annotate the graph
        xy = (new_freq,amplitude_subpix)
        xytext = (new_freq + 0.5,amplitude_subpix-0.0003)
        
        ax[1].annotate(r'$f = '+str("{:.3f}".format(new_freq))+' \,{Hz}$',xy, xytext, fontsize = text_size)
    ax[1].grid()
    
    ax[1].plot(main_freq,main_amplitude,'o', color = 'tab:blue', mec = 'k', ms = 8)
           
    plt.tight_layout()
    
    plt.savefig(fig_path + '.jpg')
    plt.savefig(fig_path + '.pdf')
    
    # write frequencies in a csv and a txt file

    suffixe = 'main_frequencies'
    sorted_indices = np.argsort(main_amplitude)
    main_freq = main_freq[sorted_indices]
    main_amplitude = main_amplitude[sorted_indices]
    
    main_peaks_propreties = (main_freq,main_amplitude)
    write_frequencies(fig_path, suffixe, signal_label, main_peaks_propreties)
    write_frequencies_csv(fig_path, suffixe, signal_label, main_peaks_propreties)
    
    return ax

def plot_multiple_features(data_dict, fps, fig_path, save, text_size = 18):
    """ Plot all features x,z,b and theta of a dictionnary data_dict and saves it. 
    
    It takes as arguments :
        - data_dict, a dictionnary of datas 
        - fps, the frame rate used during the experiment
        - fig_path, the path where the subplot will be saved
        - save, a variable to choose whether to save the plot or not, save = 1
        implies that the subplot is saved """
    
    plt.rcParams.update({'font.size': text_size})
    
    time = np.arange(0,len(data_dict['z']))/fps
    fig, ax = plt.subplots(2,2, figsize = (12,8))
    ax[0,0].plot(time,data_dict['x'] - data_dict['x'][0], color = 'tab:blue')
    ax[0,0].set_xlabel(r'Time (s)')
    ax[0,0].set_ylabel(r'$ x - x_{t=0} $ (mm)')
    ax[0,0].set_title(r'Displacement $x$ ')
    
    ax[0,1].plot(time, data_dict['z'] - np.mean(data_dict['z']), color = 'tab:orange')
    ax[0,1].set_xlabel(r'Time (s)')
    ax[0,1].set_ylabel(r'$z - \langle z \rangle$ (mm)')
    ax[0,1].set_title(r'Displacement $z$ ')
    
    ax[1,0].plot(time, data_dict['b'], color = 'tab:green')
    ax[1,0].set_xlabel(r'Time (s)')
    ax[1,0].set_ylabel(r'$b$ (mm)')
    ax[1,0].set_title(r'Minor semi-axis $b$')
    
    ax[1,1].plot(time, data_dict['theta'], color = 'tab:red')
    ax[1,1].set_xlabel(r'Time (s)')
    ax[1,1].set_ylabel(r'$\theta$ (rad)')
    ax[1,1].set_title(r'Angle $\theta$')
    
    plt.tight_layout()
    
    if save :
        plt.savefig(fig_path + '.pdf')
        plt.savefig(fig_path + '.jpg')
        write(ax, fig_path + '.pkl')
        
def plot_multiple_fourier_features(data_dict, f_sampling, fig_path, save, text_size = 18):   
    """ Plot the Fourier transform of each features x, z, b and theta and saves the subplot 
    It takes as argument :
        - data_dict, the dictionnary containing all datas 
        - f_sampling, the frequency at which the datas are sampled 
        - fig_path, the path at which the subplot is saved
        - save, a variable to choose whether the subplot is saved or not, save = 1 
        implies the subplot to be saved"""
     
    plt.rcParams.update({'font.size': text_size})
        
    fft_dict = Fourier_features(data_dict,f_sampling) 
    
    fig, ax = plt.subplots(2,2, figsize = (12,8))
    ax[0,0].plot(fft_dict['freq'],fft_dict['x'], color = 'tab:blue')
    ax[0,0].set_xlabel(r'Frequency (Hz)')
    ax[0,0].set_ylabel(r'Amplitude')
    ax[0,0].set_title(r'Displacement $x$ ')
    ax[0,0].grid()
    
    ax[0,1].plot(fft_dict['freq'], fft_dict['z'], color = 'tab:orange')
    ax[0,1].set_xlabel(r'Frequency (Hz)')
    ax[0,1].set_ylabel(r'Amplitude')
    ax[0,1].set_title(r'Displacement $z$ ')
    ax[0,1].grid()
    
    ax[1,0].plot(fft_dict['freq'], fft_dict['b'], color = 'tab:green')
    ax[1,0].set_xlabel(r'Frequency (Hz)')
    ax[1,0].set_ylabel(r'Amplitude')
    ax[1,0].set_title(r'Minor semi-axis $b$')
    ax[1,0].grid()
    
    ax[1,1].plot(fft_dict['freq'], fft_dict['theta'], color = 'tab:red')
    ax[1,1].set_xlabel(r'Frequency (Hz)')
    ax[1,1].set_ylabel(r'Amplitude')
    ax[1,1].set_title(r'Angle $\theta$')
    ax[1,1].grid()
    
    plt.tight_layout()
    
    if save :
        plt.savefig(fig_path + '.pdf')
        plt.savefig(fig_path + '.jpg')
        write(ax, fig_path + '.pkl')


def detect_peaks_features(fft_dict,min_height,window_frequency):
    
    """ Computes the main peak amplitude of all FFT spectrum of a given dictionnary. 
    It takes as arguments : 
        - fft_dict, the dictionnary containing the FFT of all features
        - min_height, the minimal height for a peak to be detected, common for all features 
        - freq_window, the frequency window within which we look for main peaks
    It returns a dictionnary containing for each features the detected peaks with their subpixel amplitude and frequency. 
    For each features, we have an array os size (2,n) n = number of detected peaks. 
    First row corresponds to frequencies and second row corresponds to their amplitude """
    
    peaks_dict = {}
    freq = fft_dict['freq'] 
    (min_freq,max_freq) = window_frequency
    
    indices_window = np.where(np.logical_and(freq > min_freq, freq < max_freq))
    freq = freq[indices_window]
    
    for key in fft_dict.keys():
        if key != 'freq':
            fourier = fft_dict[key]
            fourier = fourier[indices_window]
            ### Find peaks of the Fourier transform
            # find peaks higher than h, separated by more than distance samples from each other
            peaks_index, peaks_propreties = signal.find_peaks(fourier, min_height, distance = 1)
       
            main_freq = np.zeros(len(peaks_index))
            main_amplitude = np.zeros(len(peaks_index))
        
            for i in range(len(peaks_index)):
                # subpix precision
                idx = peaks_index[i]
                p = fourier[idx]
                p_right = fourier[idx+1]
                p_left = fourier[idx-1]
            
                delta = (p_right - p_left)/(2*(2*p - p_right - p_left))
                new_freq = freq[idx] + (freq[idx+1] - freq[idx])*delta
                
                amplitude_subpix = p + 0.5*(2*p - p_right - p_left)*delta**2
                # keep the main frequencies and amplitude
                main_freq[i] = new_freq
                main_amplitude[i] = amplitude_subpix
        
            # Sort frequencies and store them in an array 
            sorted_indices = np.argsort(main_amplitude)[::-1]
            main_freq = main_freq[sorted_indices]
            main_amplitude = main_amplitude[sorted_indices]
        
            peaks_dict[key] = np.asarray((main_freq,main_amplitude))

    return peaks_dict
    

def dephasage(s1,s2,window_frequency,f_sampling,padding_length = 2**15):
    
    """ Computes the dephasage between two signals s1 and s2. 
    The function takes as arguments :
    - the two signals of same size s1 and s2
    - window_frequency, a tuple containing values within which we look for the main peak of the FFT
    - f_sampling, the frequency at which the signals have been sampled """
    
    s1 = s1 - np.mean(s1)
    s2 = s2 - np.mean(s2)

    N = len(s1)
    # zero padding 
    if padding_length < N :
        power_2 = int(np.log(N)/np.log(2))
        padding_length = 2**(power_2 + 3)
        print('Padding length is too small for this signal')

    fourier_s1 = np.fft.fft(s1, n = padding_length, norm = 'forward')
    fourier_s2 = np.fft.fft(s2, n = padding_length, norm = 'forward')

    # computes dephasage between the two signals
    product = fourier_s1 * np.conjugate(fourier_s2)
    phase = np.angle(product[:padding_length//2], deg = False)
    # computes associated frequencies
    freq = np.fft.fftfreq(padding_length, d = 1/f_sampling)
    # keep only positive frequencies 
    freq_pos = freq[:padding_length//2]

    #keep only frequencies around the resonance frequency 
    (f_min,f_max) = window_frequency

    relevant_indices = np.where(np.logical_and(freq_pos > f_min,freq_pos < f_max))
    f_plot = freq_pos[relevant_indices]
    phase_plot = phase[relevant_indices]

    # detect main peak on s1(t)
    amplitude = np.abs(fourier_s1[:padding_length//2])
    amplitude = amplitude[relevant_indices]

    peaks_index, peaks_propreties = signal.find_peaks(amplitude, height = 0.00001, distance = 1)

    # peaks_index in descending order
    descending_indices = np.argsort(peaks_propreties['peak_heights'])[::-1]
    peaks_index = peaks_index[descending_indices]
    main_peak_idx = peaks_index[0]

    # get the phase associated to the main peak
    main_phase = phase_plot[main_peak_idx]
    main_freq = f_plot[main_peak_idx]
    
    return (main_freq,main_phase)

def compute_phase(s, window_freq, f_sampling, padding_length = 2**15):
    
    """ Computes the phase of a signal s for a frequency peak within a given window of frequencies. 
    It returns a tuple (main_frequency , main_phase) and takes as arguments : 
        - the signal s
        - window_freq, the window frequency within which we look for a peak on the Fourier spectrum 
        - f_sampling, the frequency at which the signal has been sampled
        - padding_length, a length to padd the signal in order to get better precision when computing the FFT """
    
    s = s - np.mean(s)
    N = len(s)
    # zero padding 
    if padding_length < N :
        power_2 = int(np.log(N)/np.log(2))
        padding_length = 2**(power_2 + 3)
        print('Padding length is too small for this signal')

    fourier = np.fft.fft(s, n = padding_length, norm = 'forward')

    # computes phase of the FFT 
    phase = np.angle(fourier[:padding_length//2], deg = False)
    # computes associated frequencies
    freq = np.fft.fftfreq(padding_length, d = 1/f_sampling)
    # keep only positive frequencies 
    freq_pos = freq[:padding_length//2]

    #keep only frequencies around the resonance frequency 
    (f_min,f_max) = window_freq

    relevant_indices = np.where(np.logical_and(freq_pos > f_min,freq_pos < f_max))
    f_plot = freq_pos[relevant_indices]
    phase_plot = phase[relevant_indices]

    # detect main peak on s1(t)
    amplitude = np.abs(fourier[:padding_length//2])
    amplitude = amplitude[relevant_indices]

    peaks_index, peaks_propreties = signal.find_peaks(amplitude, height = 0.00001, distance = 1)

    # peaks_index in descending order
    descending_indices = np.argsort(peaks_propreties['peak_heights'])[::-1]
    peaks_index = peaks_index[descending_indices]
    main_peak_idx = peaks_index[0]

    # get the phase associated to the main peak
    main_phase = phase_plot[main_peak_idx]
    main_freq = f_plot[main_peak_idx]
    
    return (main_freq,main_phase)