# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 08:02:28 2023

@author: sebas
"""
# Sebastien Kuchly 
#
# Definition of the different functions useful to process a video and track
# blue caps 
#
#

### Importation of useful python modules
import matplotlib.pyplot as plt
import scipy.ndimage.measurements as meas
import glob
import os
import numpy as np
import cv2 as cv
import scipy.signal as signal
from skimage.measure import label, regionprops
from scipy import ndimage
import pandas as pd 
import trackpy as tp
from skimage import util


#%% 

##################### FUNCTIONS DEFINITION ##############################

def create_dictionnary(video_file,nb_particles,fps,step):
    """ Create a dictionnary with following keywords :
        - 'time' : an array of time (in second)
        - 'particles_position' : a 3D array, dimension 1: time, dimension 2: parameters (xc,yc), 
        dimension 3: nb_particle """
        
    video_path = glob.glob(video_file)
    dt = 1/fps
    nb_frames = int(len(video_path)/step)
    particle_array = np.zeros((nb_frames,2,nb_particles))
    tracking_dict = {'time' : np.linspace(0,len(video_path)*dt,nb_frames), 
                     'particles_position' : particle_array}
    
    return tracking_dict


def particle_tracking(video_file,nb_particles,min_area,lower_blue,upper_blue):
    """ Process a video file, using trackpy module to follow several particles
    - video_file : path of the image sequence
    - nb_particles : number of particles to be detected
    - min_area : minimum pixel area of one particle 
    - lower_blue : numpy array containing the minimum blue intensity to be set to 255
    - upper_blue : numpy array containing the maximum blue intensity to be set to 255 
    This function returns a pandaDataFrame containing the following headers : 
        'x' : position along the horizontal 
        'y' : position along the vertical 
        'frame' : frame at which the position was computed 
        'particle' : index of each particles, starting from 0 """
        
    video_path = glob.glob(video_file)
    video_path.sort()
    
    # binarization of all frames
    frames = preprocessing(video_path,lower_blue,upper_blue)
    # initial DataFrame
    features = pd.DataFrame(columns = ['x', 'y', 'frame'])
    (NY,NX) = frame_dimensions(video_file)
    # loop over all frames
    for num, img in enumerate(frames):
        label_image = label(img) # we can add  background
        # intermediate DataFrame for the current frame
        intermediate_features = pd.DataFrame(columns = ['x','y','frame'])
        particles_counter = 0
        # loop over all objects
        for region in regionprops(label_image):
            # Keep only biggest particles
            if region.area > min_area :
                particles_counter += 1
                new_features = pd.DataFrame([{'y': region.centroid[0],
                                         'x': region.centroid[1],
                                         'frame': num,
                                         }])
                intermediate_features = pd.concat([intermediate_features,new_features], ignore_index = True)
        
        # particular case if some particles are not separated
        if particles_counter < nb_particles :
            # erode the frame
            new_img = frame_erosion(img, nb_particles)
            # new detection sequence with the newly eroded frame 
            label_image = label(new_img)
            intermediate_features = pd.DataFrame(columns = ['x','y','frame'])
            
            for region in regionprops(label_image):
            # Keep only biggest particles
                if region.area > min_area :
                    particles_counter += 1
                    new_features = pd.DataFrame([{'y': region.centroid[0],
                                         'x': region.centroid[1],
                                         'frame': num,
                                         }])
                    intermediate_features = pd.concat([intermediate_features,new_features], ignore_index = True)
            
        features = pd.concat([features,intermediate_features], ignore_index = True)

    # follow all trajectories
    search_range = 20
    track_features = tp.link_df(features, search_range, memory=0)
    track_features.x = track_features.x
    track_features.y = track_features.y
    
    return track_features


def DataFrame_to_dictionnary(data_frame,fps):
    """ Transforms a pandas DataFrame to a dictionnary"""
    
    time = np.array(data_frame.frame[data_frame.particle == 0]/fps)
    nb_particles = np.max(data_frame.particle) + 1
    particles_array = np.zeros((len(time),2,nb_particles))
    
    for i in range(nb_particles):
        particles_array[:,0,i] = data_frame.x[data_frame.particle == i]
        particles_array[:,1,i] = data_frame.y[data_frame.particle == i]
    
    tracking_dict = {'time' : time, 
                     'particles_position' : particles_array}
    
    return tracking_dict
    

def process_video(video_file,nb_particles,fps): 
    """ Process a video file, using the following arguments :
        - path of the video file
        - number of particles 
        - fps of the video 
        
        It returns a dictionnary as built in the function create_dictionnary """
        
    my_dict = create_dictionnary(video_file, nb_particles, fps, 1)
    
    video_path = glob.glob(video_file)
    video_path.sort()
    
    lower_blue = np.array([110,50,50])
    upper_blue = np.array([130,255,255])
    min_area = 1000
    
    frame_idx = 0
   
    for picture in video_path : 
        frame = cv.imread(picture)
    
        binarized_frame = blue_threshold(frame,lower_blue,upper_blue)
        
        particles_array = my_dict["particles_position"][frame_idx,:,:]
        my_dict["particles_position"][frame_idx,:,:] = detect_particle(binarized_frame, min_area, particles_array)
        
        frame_idx += 1
    
    return my_dict

def process_video_specified_step(video_file,nb_particles,fps,step): 
    """ Process a video file, using the following arguments :
        - path of the video file
        - number of particles 
        - fps of the video 
        - step, the time step with which each position is computed
        
        It returns a dictionnary as built in the function create_dictionnary """
        
    my_dict = create_dictionnary(video_file, nb_particles, fps, step)
    video_path = glob.glob(video_file)
    video_path.sort()
    # frames_list = get_frames_list(video_file)
    
    lower_blue = np.array([110,50,50])
    upper_blue = np.array([140,255,255])
    min_area = 2000
    
    frame_idx = 0
   
    for i in range(int(len(video_path)/step)) : 
        frame = cv.imread(video_path[i*step])
        # print(i*step)
        binarized_frame = blue_threshold(frame,lower_blue,upper_blue)
        
        particles_array = my_dict["particles_position"][frame_idx,:,:]
        my_dict["particles_position"][frame_idx,:,:] = detect_particle(binarized_frame, min_area, particles_array)
        
        frame_idx += 1
    
    return my_dict

def get_frames_list(video_file):
    """ Get the path of the video_file, format of the images should be 
        correctly precised in the video_file name such as : 
            video_file = "Z:/skuchly/Pictures/*.tiff"
            
        Then places all images in a list, which is returned by the function
    """
    video_path = glob.glob(video_file)
    video_path.sort()
    frames_list = []
    for picture in video_path :
        frame = cv.imread(picture)
        frames_list.append(frame)
    return frames_list


def blue_threshold(frame,lower_blue,upper_blue):
    """ Binarize the image with respect to the blue colour, lower_blue and 
    upper_blue corresponds to the limiting values within which a pixel will be 
    set to 255, 0 otherwise. lower_blue and upper_blue should be 1D numpy arrays in 
    the BGR system"""
    
    im_hsv = cv.cvtColor(frame, cv.COLOR_BGR2HSV)
# Threshold the HSV image to get only blue colors
    binarized_img = cv.inRange(im_hsv, lower_blue, upper_blue) 
    # filtered_image = ndimage.median_filter(binarized_img, size = 4)
    return binarized_img

def preprocessing(video_path,lower_blue,upper_blue):
    """ Binarize and erode a list of frames, it returns a list of modified frames """
    frames = []
    kernel = np.ones((5,5),np.uint8)
    for img in video_path:
        frame = cv.imread(img)
        binarized_frame = blue_threshold(frame,lower_blue,upper_blue)
        erode_img = cv.erode(binarized_frame, kernel, iterations = 1)
        frames.append(util.img_as_int(erode_img))
    return frames

def frame_erosion(frame,nb_particles):
    """ Erode a frame until all particles are detected, it returns a frame with the correctly eroded frame """
    
    particle_idx = 0
    kernel = np.ones((3,3),np.uint8)
    erosion_iteration = 0
    
    while particle_idx < nb_particles : 
        
        particle_idx = 0
        new_frame = cv.erode(frame, kernel, iterations = erosion_iteration)
        label_image = label(new_frame)
        
        # Counts the number of particles that can be observed on the frame
        for region in regionprops(label_image):
            if region.area > min_area:
                particle_idx += 1
                
        erosion_iteration += 1
    
        if erosion_iteration > 30 :
            print("Can't find particles")
            break
            
    return new_frame

def detect_particle(frame,min_area,particles_array): 
    """ Detects particles on a binarized frame. 
    Takes as argument a 2D array (nb_parameters x nb_particles)
    each particle area should be larger than the input min_area
    
    It fills the particles_array with the position of each particle"""
    
    
    kernel = np.ones((3,3),np.uint8) # array used to erode the frame
    erosion_iteration = 1
    particle_idx = 0
    nb_particles = np.shape(particles_array)[1]
    
    ## Verify if all the particles can be detected 
    while particle_idx < nb_particles : 
        
        particle_idx = 0
        new_frame = cv.erode(frame, kernel, iterations = erosion_iteration)
        label_image = label(new_frame)
        
        # Counts the number of particles that can be observed on the frame
        for region in regionprops(label_image):
            if region.area > min_area:
                particle_idx += 1
                
        erosion_iteration += 1
    
        if erosion_iteration > 30 :
            print("Can't find particles")
            break
    
    # if erosion_iteration != 2 :
    #     print("Not all particles have been detected, image needed to be eroded")
        
    ## Particle detection once the image is sufficiently eroded
    particle_idx = 0
    
    for region in regionprops(label_image):
        if region.area > min_area:
            
            circle = region
            xmin = np.min(circle.coords[:,0])
            xmax = np.max(circle.coords[:,0])
            ymin = np.min(circle.coords[:,1])
            ymax = np.max(circle.coords[:,1])
            
            # Image containing only the particle
            # We use the non-eroded image
            img_particle = frame[xmin:xmax,ymin:ymax]
            
            # Subpixel position of the particle
            (x_center, y_center) = ndimage.center_of_mass(img_particle)
            x_subpix = x_center + xmin
            y_subpix = y_center + ymin
            
            particles_array[0,particle_idx] = x_subpix
            particles_array[1,particle_idx] = y_subpix
            
            particle_idx += 1
        
    return particles_array

def frame_dimensions(video_file):
    video_path = glob.glob(video_file)
    frame = cv.imread(video_path[0])
    (NX,NY) = np.shape(frame[:,:,0])
    return (NX,NY)

def plot_displacement(tracking_dict, videofile, step, save_file):
    """ Plot the evolution of the first particle x-position with time
    x = horizontal axis
    y = vertical axis
    
    """
    
    time = tracking_dict['time']
    nb_particles = np.shape(tracking_dict['particles_position'])[2]
    (NY,NX) = frame_dimensions(videofile)
    frequency = tracking_dict['frequency'] 
    
    # fps used to analyze 
    fps_tracking = 1+int(1/(time[1] - time[0]))
    
    for i in range(nb_particles):
    
        xc = tracking_dict['particles_position'][:,0,i]
        yc = tracking_dict['particles_position'][:,1,i]

        fig, ax = plt.subplots(2,1)
        fig.suptitle('Position evolution of particle #'+str(i)+', fps = '+str(fps_tracking)+', frequency = '+str("{:.2f}".format(frequency))+'Hz')
        fig.subplots_adjust(wspace = 0, hspace = 0.7)
        
        ax[0].plot(time,xc)
        ax[0].set_xlabel('t (s)')
        ax[0].set_ylabel('x')
        ax[0].set_title('x(t)')
        ax[0].grid()
        
        ax[1].plot(time,yc)
        ax[1].set_xlabel('t (s)')
        ax[1].set_ylabel('y')
        ax[1].set_title('y(t)')
        ax[1].grid()
    
        plt.savefig(save_file + 'Position_evolution_particle_n' + str(i)+ '_fps' +str(fps_tracking)+ '_frequency'+str("{:.2f}".format(frequency))+'Hz.png')
    
def get_velocity(tracking_dict):
    """ Computes the velocity (in pixels/s) of each particles (at first order) 
    It creates within the given dictionnary, new keys :
        - 'particles_velocity' : a 3D array, dimension 1: time, dimension 2: parameters (u,v), 
        dimension 3: nb_particle """
    
    time = tracking_dict['time']
    dt = time[1] - time[0]
    nb_particles = np.shape(tracking_dict['particles_position'])[2]
    velocity_array = np.zeros((len(time)-1,2,nb_particles))
    
    
    for i in range(nb_particles):
        x = tracking_dict['particles_position'][:,0,i]
        y = tracking_dict['particles_position'][:,1,i]
        velocity_array[:,0,i] = (x[1:] - x[:-1])/dt
        velocity_array[:,1,i] = (y[1:] - y[:-1])/dt
    
    
    tracking_dict['particles_velocity'] = velocity_array
    # tracking_dict['velocity_time'] = time[:-1]
    
    # velocity_dict = {'time' : time[:-1], 
    #                  'particles_velocity' : velocity_array}
    return tracking_dict

    
def plot_velocity(tracking_dict, scale, step, save_file):
    """ Plot the evolution of the particles velocity field with time
    The velocity field is expressed in cm/s  
    
    1 pixel = 1 scale cm
    """
    
    time = tracking_dict['time'][:-1]
    nb_particles = np.shape(tracking_dict['particles_velocity'])[2]
    frequency = tracking_dict['frequency']
    # fps used to analyze 
    fps_tracking = 1+int(1/(time[1] - time[0]))
    
    for i in range(nb_particles):
    
        u = tracking_dict['particles_velocity'][:,0,i]*scale
        v = tracking_dict['particles_velocity'][:,1,i]*scale
    
        fig, ax = plt.subplots(2,1)
        fig.suptitle('Velocity field of particle #'+str(i)+', fps = '+str(fps_tracking)+', frequency = '+str("{:.2f}".format(frequency))+'Hz')
        fig.subplots_adjust(wspace = 0, hspace = 0.7)
        
        ax[0].plot(time,u)
        ax[0].set_xlabel('t (s)')
        ax[0].set_ylabel('u')
        ax[0].set_title('u(t) (cm/s)')
        ax[0].grid()
        
        ax[1].plot(time,v)
        ax[1].set_xlabel('t (s)')
        ax[1].set_ylabel('v')
        ax[1].set_title('v(t) (cm/s)')
        ax[1].grid()
        
        plt.savefig(save_file + 'Velocity_evolution_particle_n' + str(i)+ '_fps' +str(fps_tracking)+ '_frequency'+str("{:.2f}".format(frequency))+'Hz.png')


def plot_2D_trajectory_xy(tracking_dict, video_file, step):
    """ We should Try to conserve particles label at each frames 
    - x corresponds to the vertical axis
    - y correspons to the horizontal axis """
    (NX,NY) = frame_dimensions(video_file)

    nb_particles = np.shape(tracking_dict['particles_position'])[2]
    fig, ax = plt.subplots()
    fig.suptitle('2D position tracking of particles #, PIV step ='+str(step)+'')
    
    colour = tracking_dict['time']
    
    for i in range(nb_particles):
        
        x = tracking_dict['particles_position'][:,0,i]/NX
        y = tracking_dict['particles_position'][:,1,i]/NY
        
        plot = ax.scatter(y,x, c = colour, cmap = 'viridis', s = 5)
        
        ax.set_ylabel('x')
        ax.set_xlabel('y')
        ax.invert_yaxis()
        
    fig.colorbar(plot)

def plot_2D_trajectory_yx(tracking_dict, video_file, step, save_file):
    """ We should Try to conserve particles label at each frames 
    - x corresponds to the horizontal axis
    - y correspons to the vertical axis """
    (NY,NX) = frame_dimensions(video_file)

    nb_particles = np.shape(tracking_dict['particles_position'])[2]
    fig, ax = plt.subplots()
    fig.suptitle('2D position tracking of particles, PIV step ='+str(step)+'')
    
    colour = tracking_dict['time']
    
    for i in range(nb_particles):
        
        x = tracking_dict['particles_position'][:,0,i]/NX
        y = tracking_dict['particles_position'][:,1,i]/NY
        
        plot = ax.scatter(x,y, c = colour, cmap = 'viridis', s = 5)
        
        ax.set_ylabel('y')
        ax.set_xlabel('x')
        ax.invert_yaxis()
        
    fig.colorbar(plot)
    
    plt.savefig(save_file + '2D_particle_tracking.png')

def Fourier_spectrum(s,f_sampling):
    """ Computes and plot the Fourier spectrum of a signal s, using the following arguments :
        - signal s
        - f_sampling = fps/step, the frequency at which the signal is sampled """
    
    N = len(s)
    s = s[:N] - np.mean(s[:N])
    fourier = np.fft.fft(s)

    # Get frequency array
    # first argument corresponds to the size of the window
    # second argument d corresponds to the sample spacing (inverse of the sampling rate)
    freq = np.fft.fftfreq(N,d = 1/f_sampling)
    
    # Take the absolute value of the Fourier transformation, only for positive frequencies
    fourier_abs = np.abs(fourier[:N//2])
    # Normalization of the amplitudes
    fourier_norm = fourier_abs*2/N
    # Keep only positive frequencies
    freq_pos = freq[:N//2]

    fig,ax = plt.subplots()
    fig.suptitle('Fourier spectrum')

    ax.plot(freq_pos,fourier_norm, label = 'Normalized amplitudes')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Amplitude')


def highpass_filter(s,signal_time,f_cut,filter_order,f_sampling):
    """ Apply a highpass filter to a signal s, using signal.butter function
    - s = the signal
    - signal_time = time serie of the signal
    - f_cut = cutting frequency 
    - filter_order = order of the filter applied 
    - f_sampling = frequency at which the signal is sampled 
    
    It returns a numpy array of the filtered signal """
    
    sos = signal.butter(filter_order,f_cut,'highpass', fs = f_sampling, output = 'sos')
    # filtered signal 
    filtered = signal.sosfilt(sos, s)
    
    # Plot of the new signal
    fig, ax = plt.subplots(2,1)
    fig.subplots_adjust(wspace = 0, hspace = 0.6)
    
    ax[0].plot(signal_time,s)
    ax[0].set_xlabel('time (s)')
    ax[0].set_ylabel('s(t)')
    ax[0].set_title('Initial signal')
    ax[0].grid()
    
    ax[1].plot(signal_time,filtered)
    ax[1].set_xlabel('time (s)')
    ax[1].set_ylabel('s(t) filtered')
    ax[1].set_title('Transformed signal')
    ax[1].grid()
    
    return filtered
    
#%%
################### MAIN PROGRAM  USING TRACKPY #######################

lower_blue = np.array([110,50,50])
upper_blue = np.array([140,255,255])
min_area = 1000
nb_particles = 13
fps = 100

path = 'W:/Banquise/Sebastien/Experiences/Waves_field_comparison_test_25_04_23/25_04_23_v100_T150_17particles/'

video_file  = path  + 'image_sequence/*.tiff'
save_file = path + 'results/'

# creation of a directory if no directory for the results
if os.path.isdir(save_file) == False:
    os.mkdir(save_file)

track_data = particle_tracking(video_file, nb_particles, min_area, lower_blue, upper_blue)

tracking_dict = DataFrame_to_dictionnary(track_data, fps)


#%%
step = 1
plot_2D_trajectory_yx(tracking_dict, video_file, step, save_file)

#%% 
tracking_dict['frequency'] = 1/300
plot_displacement(tracking_dict, video_file, step, save_file)

#%% 
##################### MAIN PROGRAM ##########################

lower_blue = [110,50,50]
upper_blue = [140,255,255]
min_area = 1000

filename = 'W:\Banquise\Sebastien\Experiences\Test_particle_tracking_videos/12_04_2023_v150_T150/*.tiff'
nb_particles = 1
fps = 100
step = 5
T_change = 150 # period at which the engine rotation is modified (in ms)

# frame_list = get_frames_list(filename)
# print(len(frame_list))

my_dict = process_video_specified_step(filename, nb_particles, fps, step)
# print(my_dict['particles_position'][:,:,0])

# plot_positions(my_dict)

time = my_dict['time']
nb_particles = np.shape(my_dict['particles_position'])[2]


#%% 

video_path = glob.glob(filename)
img = cv.imread(video_path[0])

plt.figure()
plt.imshow(img)
#%%

### Plot of x(t) and y(t)

T_change = 150 # period at which the engine rotation is modified (in ms)
my_dict['frequency'] = 1000/(2*T_change) # frequency used in the experiment (in Hz)
plot_displacement(my_dict, filename,step)

#%% 

plot_2D_trajectory(my_dict, filename)

#%%

scale = 0.04 # 1 pixel = scale cm

velocity_dict = get_velocity(my_dict)
velocity_dict['frequency'] = 1000/(2*T_change)
plot_velocity(velocity_dict,scale,step)

#%%

x = velocity_dict['particles_velocity'][:,0,0]
signal_time = velocity_dict['time'][:-1]
f_cut = 0.5
filter_order = 10
f_sampling = fps/step
x_filter = highpass_filter(x, signal_time, f_cut, filter_order, f_sampling)
Fourier_spectrum(x, f_sampling)
Fourier_spectrum(x_filter, f_sampling)

#%%

############## Try to get the Fourier components of the particle trajectory #########
#### Be careful to take the correct coordinate, parallel to waves

(NX,NY) = frame_dimensions(filename)

x = my_dict['particles_position'][:,0,0]/NX
y = my_dict['particles_position'][:,1,0]/NY
time = my_dict['time']

N = len(x)
s = x[:N] - np.mean(x[:N])
signal_time = time[:N]
fourier = np.fft.fft(s)

# test = np.sin(time[:200])
# fourier = np.fft.fft(test)

fig, ax = plt.subplots(3,1)

ax[0].plot(signal_time,s,'.',ms = 1)
ax[1].plot(fourier.real)
ax[2].plot(fourier.imag)

# Get frequency array
# first argument corresponds to the size of the window
# second argument d corresponds to the sample spacing (inverse of the sampling rate)
freq = np.fft.fftfreq(len(signal_time),d = step/fps)

fig, ax = plt.subplots(2,1)

ax[0].plot(freq, fourier.real, label = 'Real part')
ax[0].plot(freq, fourier.imag, label = 'Imaginary part')
ax[0].legend()

ax[1].plot(freq, fourier.real, label = 'Real part')
ax[1].plot(freq, fourier.imag, label = 'Imaginary part')
ax[1].set_xlim([-2.5,2.5])
ax[1].legend()


#%%
#### We try to keep only positive frequencies and normalized amplitude 

# Take the absolute value of the Fourier transformation, only for positive frequencies
fourier_abs = np.abs(fourier[:N//2])
# Normalization of the amplitudes
fourier_norm = fourier_abs*2/N
# Keep only positive frequencies
freq_pos = freq[:N//2]

fig,ax = plt.subplots()

ax.plot(freq_pos,fourier_norm, label = 'Normalized amplitudes')
ax.set_xlabel('Frequency (Hz)')
ax.set_ylabel('Amplitude')


#%%
#### Spectrogramm of a measurement x(t)

# Calcul du spectrogramme
f, t, Sxx = signal.spectrogram(s, fs = fps/step)

# On limite aux fréquences présentent
Sxx_red = Sxx[np.where(f<100)]
f_red = f[np.where(f<100)]

# Affichage du spectrogramme
plt.figure()
plt.pcolormesh(t, f_red, Sxx_red, shading = 'gouraud')
plt.ylabel('Frequency (Hz)')
plt.xlabel('Time (s)')
plt.title('Spectrogramm')
plt.show()

#%% 

#### Fourier transformation using scipy.signal.butter

(NX,NY) = frame_dimensions(filename)

x = my_dict['particles_position'][:,0,0]/NX
y = my_dict['particles_position'][:,1,0]/NY
time = my_dict['time']

N = len(x)
s = x[:N] - np.mean(x[:N])
signal_time = time[:N]
fourier = np.fft.fft(s)

# 1st plot
fig, ax = plt.subplots(3,1)

ax[0].plot(signal_time,s,'.',ms = 1)
ax[1].plot(fourier.real)
ax[2].plot(fourier.imag)

# Get frequency array
# first argument corresponds to the size of the window
# second argument d corresponds to the sample spacing (inverse of the sampling rate)
freq = np.fft.fftfreq(len(signal_time),d = step/fps)

# 2nd plot
fig, ax = plt.subplots(2,1)

ax[0].plot(freq, fourier.real, label = 'Real part')
ax[0].plot(freq, fourier.imag, label = 'Imaginary part')
ax[0].legend()

ax[1].plot(freq, fourier.real, label = 'Real part')
ax[1].plot(freq, fourier.imag, label = 'Imaginary part')
ax[1].set_xlim([-2.5,2.5])
ax[1].legend()

# apply a highpass to the signal

filter_order = 10
freq_cut = 0.5
f_sampling = fps/step
sos = signal.butter(filter_order,freq_cut,'highpass', fs = f_sampling, output = 'sos')
filtered = signal.sosfilt(sos, s)

fig, ax = plt.subplots(2,1)

ax[0].plot(signal_time,s)
ax[1].plot(signal_time,filtered)


