# -*- coding: utf-8 -*-
"""
Created on Mon Dec 23 10:01:39 2024

@author: sebas
"""



import math
import numpy as np
from scipy.fftpack import fft, ifft


#----------------------------------------------------------------------------------------------------

def nextpow2(x):
    """ Return closer exposant p such that 2**p > x
    """
    if x == 0:
        p = 0
    elif x == 1:
        p = 1
    else:
        p = math.ceil(math.log2(x))
    return p

#---------------------------------------------------------------------------------------------------

def fft_1D(s,N,fs):
    """ Compute FFT of a 1D signal 
    Inputs : 
        - s : 1D array,
        - N : size of the array after padding,
        - fs : sampling frequency """
        
    original_length = len(s)
    FFT = fft(s - np.mean(s),n = N)
    FFT = FFT[:N//2]/original_length
    FFT[1:-1] = 2*FFT[1:-1]
    
    freq = fs*np.arange(0,(N//2))/N
    
    return FFT,freq


#---------------------------------------------------------------------------------------------------

def temporal_FFT(H,fps,padding_bool = 1,add_pow2 = 0,output_FFT = False):
    
    original_length = np.size(H,2) # original length of the signal
    padding_length = 2**(nextpow2(original_length) + add_pow2)

    # averaging over time 
    H_mean = np.mean(H,2)
    H_mean = np.tile(H_mean,(np.size(H,2),1,1))
    H_mean = np.transpose(H_mean,(1,2,0))
    H = H - H_mean

    if padding_bool :
        FFT_t = np.fft.fft(H,n = padding_length,axis = 2)
        N = padding_length 
        print('Padding used')
    else :
        FFT_t = np.fft.fft(H,n = N,axis = 2)
        N = original_length
        print('No padding')


    FFT_t = FFT_t/original_length # normalization of the FFT
    TF_inter = np.mean(abs(FFT_t),(0,1))

    TF_spectrum = TF_inter[:N//2]
    TF_spectrum[1:-1] = 2*TF_spectrum[1:-1] # multiply by 2 for peaks that are both in positive an negative frequencies
    # freq = np.fft.fftfreq(N,d = 1/fps)
    freq = fps*np.arange(0,(N/2))/N
    # freq = freq[:N//2]
    
    if not output_FFT :
        return TF_spectrum,freq
    
    else : 
    # keep only one half of the spectrum
        FFT_t = FFT_t[:,:,:N//2]
        FFT_t[:,:,1:-1] = 2*FFT_t[:,:,1:-1]
    
        return TF_spectrum,freq,FFT_t

#----------------------------------------------------------------------------------------------------

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


#--------------------------------------------------------------------------------------------------------

def polyfit2d(x, y, z, n=3, m=3, order=None):
    '''
    Two dimensional polynomial fitting by least squares.
    Fits the functional form f(x,y) = z.

    Notes
    -----
    Resultant fit can be plotted with:
    np.polynomial.polynomial.polygrid2d(x, y, soln.reshape((n+1, m+1)))

    Parameters
    ----------
    x, y: array-like, 1d
        x and y coordinates.
    z: np.ndarray, 2d
        Surface to fit.
    n, m: int, default is 3
        Polynomial order in x and y, respectively.
    order: int or None, default is None
        If None, all coefficients up to maxiumum n, m, ie. up to and including x^n*y^m, are considered.
        If int, coefficients up to a maximum of n+m <= order are considered.

    Returns
    -------
    Return paramters from np.linalg.lstsq.

    soln: np.ndarray
        Array of polynomial coefficients.
    residuals: np.ndarray
    rank: int
    s: np.ndarray

    '''

    # grid coords
    x, y = np.meshgrid(x, y)
    # coefficient array, up to x^n, y^m
    coeffs = np.ones((n+1, m+1))

    # solve array
    a = np.zeros((coeffs.size, x.size))

    # for each coefficient produce array x^i, y^j
    for index, (i, j) in enumerate(np.ndindex(coeffs.shape)):
        # do not include powers greater than order
        if order is not None and i + j > order:
            arr = np.zeros_like(x)
        else:
            arr = coeffs[i, j] * x**i * y**j
        a[index] = arr.ravel()

    # do leastsq fitting and return leastsq result
    soln, residuals, rank, s = np.linalg.lstsq(a.T, np.ravel(z), rcond=None)
    return soln, residuals, rank, s

#--------------------------------------------------------------------------------------------------------------

def supress_quadratic_noise(V,x,y):
    """ Supress noise associated to drone motion. Drone motion is supposed to introduce a quadratic field to 
    the wave field measured using DIC. 
    Inputs : - V, velocity field, np.ndarray 3D [nx,ny,nt]
             - x, x-coordinates, np.ndarray 1D
             - y, y-coordinates, np.ndarray 1D
             
    Outputs : - Vs, velocity field corrected from quadratic noise, np.ndarray 3D [nx,ny,nt] """
    
    Vs = np.zeros(np.shape(V))
    nx = 2
    ny = 2
    
    X,Y = np.meshgrid(x,y,copy = False)
    
    print('Supressing quadratic noise...')
    for k in range(np.shape(V)[2]):
        field = V[:,:,k]
        # Fit current field by a quadratic polynom using least square method 
        soln,residuals,rank,s = polyfit2d(x,y,field,n = nx,m = ny)
        fitted_field = np.polynomial.polynomial.polyval2d(X, Y, soln.reshape((nx+1,ny+1)))
        
        Vs[:,:,k] = field - fitted_field
    
    print('Quadratic noise supressed')
    return Vs

#------------------------------------------------------------------------------------------------------------------

def radial_average(matrix, center):
    """
    Compute the radial average of a 2D matrix with respect to a specified center.
    
    Parameters:
        matrix (2D array): Input 2D matrix (ny,nx).
        center (tuple): Coordinates (x0, y0) of the center.
    
    Returns:
        radii (1D array): Unique radii values.
        radial_profile (1D array): Radially averaged intensity values.
    """
    # Create a grid of indices
    y, x = np.indices(matrix.shape)
    x0, y0 = center

    # Compute radial distances from the center
    r = np.sqrt((x - x0)**2 + (y - y0)**2)

    # Flatten the matrix and distances for binning
    r_flat = r.flatten()
    values_flat = matrix.flatten()

    # Get sorted indices based on the radial distances
    sort_indices = np.argsort(r_flat)
    r_sorted = r_flat[sort_indices]
    values_sorted = values_flat[sort_indices]

    # Find unique radii and their indices
    unique_radii, radii_indices = np.unique(r_sorted, return_index=True)

    # Compute radial averages
    radial_profile = [
        values_sorted[start:end].mean() 
        for start, end in zip(radii_indices[:-1], radii_indices[1:])
    ]

    # Append the last average
    radial_profile.append(values_sorted[radii_indices[-1]:].mean())

    return unique_radii, np.array(radial_profile)

#---------------------------------------------------------------------------------------------------------------

def radial_average_bin(matrix, center, bin_size=1.0):
    """
    Compute the radial average of a 2D matrix with respect to a specified center using larger steps (bins).
    
    Parameters:
        matrix (2D array): Input 2D matrix (image).
        center (tuple): Coordinates (x0, y0) of the center.
        bin_size (float): Size of the radial bins for averaging.
    
    Returns:
        bin_edges (1D array): Edges of the radial bins.
        radial_profile (1D array): Radially averaged intensity values for each bin.
    """
    # Create a grid of indices
    y, x = np.indices(matrix.shape)
    x0, y0 = center

    # Compute radial distances from the center
    r = np.sqrt((x - x0)**2 + (y - y0)**2)

    # Flatten the matrix and distances for binning
    r_flat = r.flatten()
    values_flat = matrix.flatten()

    # Define bins based on the bin size
    max_radius = np.max(r_flat)
    bin_edges = np.arange(0, max_radius + bin_size, bin_size)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Digitize the radial distances into bins
    bin_indices = np.digitize(r_flat, bin_edges) - 1  # Subtract 1 to make it zero-based

    # Compute radial averages for each bin
    radial_profile = [
        values_flat[bin_indices == i].mean() if np.any(bin_indices == i) else 0 
        for i in range(len(bin_centers))
    ]

    return bin_centers, np.array(radial_profile)

#-------------------------------------------------------------------------------------------------------------

def space_time_spectrum(V,facq_x,facq_t,add_pow2 = [0,0,0]):
    """" Compute space-time spectrum E(f,k) 
    Inputs : - V, np.ndarray  dim : [ny,nx,nt]
             - facq_x, float, acquisition frequency in space 
             - facq_t, float, acquisition frequency in time
             - add_pow2, np.ndarray or list, additional padding 
             
    Outputs : a dictionary with the following keys : 
              - E, space-time spectrum, np.ndarray dim : [nf,nk] 
              - shift, space-time FFT, with only positive frequencies, np.ndarray, dim : [nky,nkx,nf]
              - k, list of wave vectors magnitude
              - freq, list of frequencies (positive)
              - kx, list of wave vectors x-coordinate
              - ky, list of wave vectors y-coordinate """

    padding = [2**(nextpow2(np.shape(V)[d]) + add_pow2[d]) for d in range(3)]
    FFT = np.fft.fftn(V,s = padding ,axes = (0,1,2))/np.size(V)
    
    # compute array of frequencies and wave vectors 
    freq = facq_t*np.arange(0,(padding[2]/2))/padding[2]
    
    kx = 2*np.pi*facq_x*np.arange(-padding[1]/2,padding[1]/2-1)/padding[1]
    ky = 2*np.pi*facq_x*np.arange(-padding[0]/2,padding[0]/2-1)/padding[0]
    
    # keep only positive frequencies 
    FFT_positive = FFT[:,:,:padding[2]//2]
    FFT_positive[:,:,1:-1] = 2*FFT_positive[:,:,1:-1]
    
    shift = np.fft.fftshift(FFT_positive,axes = (0,1))
    
    # center of the FFT 
    x0 = padding[1]//2
    y0 = padding[0]//2
    
    # radial average over wave vector directions
    E = []
    for idx in range(np.shape(shift)[2]):
        field = shift[:,:,idx]
        R_tics, R_profile = radial_average_bin(abs(field), (x0,y0))
        E.append(R_profile)
        
    k = 2*np.pi*R_tics*facq_x/padding[0]
    E = np.array(E)
    
    result = {'E':E,'shift':shift,'k':k,'freq':freq,'kx':kx,'ky':ky}
    
    return result
