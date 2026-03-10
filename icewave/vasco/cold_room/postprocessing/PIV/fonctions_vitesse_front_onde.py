# %% import libraries
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.interpolate import interp2d
from scipy.interpolate import interp1d
import h5py
import os
from skimage.measure import profile_line
from matplotlib.animation import FuncAnimation
from scipy import signal
from scipy.optimize import curve_fit

# %% fonctions pour importer les données et les convertir
def extract_data_mat_file(name_mat_file):
    with h5py.File(name_mat_file, 'r') as file:
        # List all the keys in the file
        #print("Keys: %s" % list(file.keys()))
        # Access a dataset
        data = file['data'][:]  # Accessing 'data' variable
    return data
# Now you can work with the 'data' variable as you would with any other Python variable

def tuples_to_complex(matrix):
    """
    Convert a 2D matrix containing tuples into a matrix of complex numbers.
    
    Args:
    matrix (numpy.ndarray): 2D matrix containing tuples (a, b).
    
    Returns:
    numpy.ndarray: 2D matrix where tuples are replaced with complex numbers of the form a + ib.
    """
    # Initialize an empty matrix to store complex numbers
    complex_matrix = np.empty(matrix.shape, dtype=complex)
    
    # Iterate through each element in the matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            # Get the tuple from the original matrix
            a, b = matrix[i, j]
            
            # Create a complex number from the tuple and assign it to the corresponding position
            complex_matrix[i, j] = complex(a, b)
    
    return complex_matrix

def convert_matfile_to_complex_field(name_mat_file):
    data = extract_data_mat_file(name_mat_file)
    complex_field = tuples_to_complex(data)
    return complex_field

def convert_matfile_to_image(name_mat_file):
    data = extract_data_mat_file(name_mat_file)
    data_complex = tuples_to_complex(data)
    image = np.real(data_complex)/np.abs(data_complex)
    return image

def compute_real_field(complex_field,t,freq):
    return np.real(complex_field * np.exp(-1j*2*np.pi*freq*t))/np.abs(complex_field)

# %% fonctions pour obtenir les coordonnées des points sur le cercle à partir des propriétés deu cercle
# les arguments importants à rentrer sont coord_source et radius

def Circle(x,O=(0,0),R=1):#coord de O dans l'ordre x_coord,y_coord (où l'axe y est inversé)
    yp = O[1] + np.sqrt(R**2 - (x-O[0])**2)
    ym = O[1] - np.sqrt(R**2 - (x-O[0])**2)
    return yp,ym

def coord_circle_on_image(complex_field,coord_source,radius,plot=False):
    image = np.real(complex_field)/np.abs(complex_field)
    xvals = np.arange(image.shape[1])
    x_circle = np.hstack((xvals,xvals))
    yp,ym = Circle(xvals,O=coord_source,R=radius)
    y_circle = np.hstack((yp,ym))

    #complex_field = convert_matfile_to_complex_field(name_mat_file)
    indices_to_plot = (x_circle>0)&(x_circle<complex_field.shape[1])&(y_circle>0)&(y_circle<complex_field.shape[0])

    x_circle_on_image = x_circle[indices_to_plot]
    y_circle_on_image = y_circle[indices_to_plot]

    if plot==True:
        plt.figure()
        plt.imshow(image)
        plt.plot(x_circle_on_image,y_circle_on_image,'r')
        plt.xlabel('x',fontsize=15)
        plt.ylabel('y',fontsize=15)
        #plt.ylim(-50,120)
        plt.show()
    
    return x_circle_on_image,y_circle_on_image,x_circle,y_circle

# %% fonctions pour calculer les profils du champ de hauteur du point source à un point donné sur le cercle
def field_along_profile(image,coord_start_point,coord_end_point,plot=False):
    image_with_nans = np.vstack((np.ones_like(image)*np.nan , image))

    #index_point_on_circle = 10
    #coord_start_point = coord_source
    #coord_end_point = (np.round(x_circle_on_image[index_point_on_circle]).astype(int),np.round(y_circle_on_image[index_point_on_circle]).astype(int))

    coord_start_point_with_nans = (coord_start_point[0],coord_start_point[1]+image.shape[0])
    coord_end_point_with_nans = (coord_end_point[0],coord_end_point[1]+image.shape[0])

    #print(image.shape[0])
    #print(coord_start_point,coord_end_point)
    #print(coord_start_point_with_nans,coord_end_point_with_nans)

    p = profile_line(image_with_nans,np.flip(coord_start_point_with_nans),np.flip(coord_end_point_with_nans))
    
    if plot==True:
        plt.figure()
        plt.imshow(image_with_nans,vmin=-1,vmax=1)
        plt.plot(coord_start_point_with_nans[0],coord_start_point_with_nans[1],'ob')
        plt.plot(coord_end_point_with_nans[0],coord_end_point_with_nans[1],'or')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.colorbar()
        plt.show()        
    return p

def field_along_profile_endpointoncircle(image,coord_source,x_circle_on_image,y_circle_on_image,index_point_on_circle,plot=False): # attention le probleme est que ca ne doit pas depasser la longueur du tableau x_circle_on_image
    coord_start_point = coord_source
    coord_end_point = (np.round(x_circle_on_image[index_point_on_circle]).astype(int),np.round(y_circle_on_image[index_point_on_circle]).astype(int))
    p = field_along_profile(image,coord_start_point,coord_end_point,plot=plot)
    return p

# certaines fonctions n'ont pas été écrites ici

############################################################
# fonction qui smooth une image ############################
############################################################
def sgolay2d ( z, window_size, order, derivative=None):
    """
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0
    
    if  window_size % 2 == 0:
        raise ValueError('window_size must be odd')
    
    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2
    
    # exponents of the polynomial. 
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
    # this line gives a list of two item tuple. Each tuple contains 
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]
    
    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    A = np.empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])
        
    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band ) 
    # left band
    band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z
    
    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band ) 
    
    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band ) 
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band ) 
    
    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return signal.fftconvolve(Z, -c, mode='valid')        
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return signal.fftconvolve(Z, -r, mode='valid')        
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return signal.fftconvolve(Z, -r, mode='valid'), signal.fftconvolve(Z, -c, mode='valid') 

# fonction pur interpoler un tableau 1d en utilisant le facteur par lequel on veut multiplier le nb de pts:
def interpolate_1d_profile(p,alpha): # alpha is the rescaling factor we want to apply for the interpolation
    yp = p
    xp = np.arange(len(yp))
    x = np.linspace(0,np.max(xp),len(xp)*alpha)
    y = np.interp(x,xp,yp)
    return x,y
    
# fonction pour interpoler un champ réel ("image")
def interpolate_real_field(real_field):
    x = np.arange(real_field.shape[1])
    y = np.arange(real_field.shape[0])
    f = interp2d(x,y,real_field,kind='cubic')
    xnew = np.linspace(0,len(x),int(len(x)))
    ynew = np.linspace(0,len(y),int(len(y)))
    real_field_interp = f(xnew,ynew)
    return real_field_interp
    
# fonctions qui correle des profils à temps successifs
def renormalize_profile(p):
    minp = np.min(p)
    maxp = np.max(p)
    if minp==0:
        minp=np.min(np.where(p==0,np.inf,p))
    if maxp==0:
        maxp=np.max(np.where(p==0,-np.inf,p))
    amplitude = (maxp - minp)/2
    shift = (maxp + minp)/2
    return p/amplitude - shift/amplitude
def correlate_2profiles(complex_field,time_before,time_after,index_point_on_circle,coord_source,x_circle_on_image,y_circle_on_image,freq,shift_startpoint=0,renormalize=False,smooth=True,window_size=41,order=4,plot=False):
    #image_a = interpolate_real_field(compute_real_field(complex_field,time_before,freq=freq))
    #image_b = interpolate_real_field(compute_real_field(complex_field,time_after,freq=freq))
    image_a = compute_real_field(complex_field,time_before,freq=freq)
    image_b = compute_real_field(complex_field,time_after,freq=freq)
    if smooth==True:
        window_size_new = int(window_size/2) * 2 + 1
        image_a = sgolay2d(image_a,window_size_new,order)
        image_b = sgolay2d(image_b,window_size_new,order)
    coord_source_new = coord_source
    x_circle_on_image_new = x_circle_on_image
    y_circle_on_image_new = y_circle_on_image
    p_a = field_along_profile_endpointoncircle(image_a,coord_source_new,x_circle_on_image_new,y_circle_on_image_new,index_point_on_circle,plot=plot)[shift_startpoint:]
    p_b = field_along_profile_endpointoncircle(image_b,coord_source_new,x_circle_on_image_new,y_circle_on_image_new,index_point_on_circle,plot=plot)[shift_startpoint:]
    p_a_with_zeros = np.where(np.isnan(p_a),0,p_a)
    p_b_with_zeros = np.where(np.isnan(p_b),0,p_b)
    if renormalize==True:
        p_a_with_zeros = renormalize_profile(p_a_with_zeros)
        p_b_with_zeros = renormalize_profile(p_b_with_zeros)
        # et on refait l'operation suivante à cause du shift
        p_a_with_zeros = np.where(np.isnan(p_a),0,p_a_with_zeros)
        p_b_with_zeros = np.where(np.isnan(p_b),0,p_b_with_zeros)
    p_a_on_image = np.delete(p_a_with_zeros,np.where(np.isnan(p_a))) ####### ATTENTION
    p_b_on_image = np.delete(p_b_with_zeros,np.where(np.isnan(p_b))) ####### ATTENTION
    correlation_a_b = signal.correlate(p_a_on_image,p_b_on_image,mode='valid')
    lags = signal.correlation_lags(len(p_a_on_image), len(p_b_on_image)) # voir comment sont calculés les lags
    return lags,correlation_a_b,p_a_on_image,p_b_on_image

def compute_dist_max_corr(complex_field,time_before,time_after,index_point_on_circle,coord_source,x_circle_on_image,y_circle_on_image,freq,shift_startpoint=0,plot=False,smooth=True,window_size=71,order=2):
    lags,correlation_a_b,p_a_with_zeros,p_b_with_zeros = correlate_2profiles(complex_field,time_before,time_after,index_point_on_circle,coord_source,x_circle_on_image,y_circle_on_image,freq,shift_startpoint=shift_startpoint,plot=plot,smooth=smooth,window_size=window_size,order=order)
    index_max_correlation_a_b = np.where(correlation_a_b==np.max(correlation_a_b))
    dist_max_corr = lags[index_max_correlation_a_b] # distance par unité de cases PIV
    if plot==True:
        plt.figure()
        plt.plot(p_a_with_zeros)
        plt.plot(p_b_with_zeros)
        plt.show()
        plt.figure()
        plt.plot(lags,correlation_a_b)
        print(dist_max_corr)
        plt.vlines(dist_max_corr,np.min(correlation_a_b),np.max(correlation_a_b),linewidth=0.7,linestyle='--',color='red')
        plt.show()
    return dist_max_corr,correlation_a_b,p_a_with_zeros,p_b_with_zeros
# code pour créer la video demodulée dans vitesse_front_onde_2
# fit qui utilise plusieurs dist_max_corr par rapport à l'instant 0
def linear(x,a):
    return a*x

#(complex_field,tab_times_after,index_point_on_circle,coord_source,x_circle_on_image,y_circle_on_image,alpha,f_exc,shift_startpoint=shift_startpoint,plot=False)
def correlate_full_period_with_t0(complex_field,tab_times_after,index_point_on_circle,coord_source,x_circle_on_image,y_circle_on_image,freq,shift_startpoint=0,plot=False,smooth=True):
    tab_dist_max_corr = np.zeros_like(tab_times_after)
    for i in range(len(tab_times_after)):
        dist_max_corr = compute_dist_max_corr(complex_field,0,tab_times_after[i],index_point_on_circle,coord_source,x_circle_on_image,y_circle_on_image,freq,shift_startpoint=shift_startpoint,plot=plot,smooth=smooth)[0]
        #print(dist_max_corr[0])
        tab_dist_max_corr[i] = dist_max_corr[0]
    return tab_dist_max_corr

def compute_phase_velocity(tab_times_after,tab_dist_max_corr,plot=False,unwrap=True):
    if unwrap==True:
        tab_dist_max_corr = np.unwrap(tab_dist_max_corr)
    #dt = tab_times_after[1] - tab_times_after[0]
    opt,pcov = curve_fit(linear,tab_times_after,tab_dist_max_corr)
    if plot==True:
        plt.figure()
        plt.plot(tab_times_after,tab_dist_max_corr,'+')
        plt.plot(tab_times_after,linear(tab_times_after,opt[0]))
        plt.xlabel('delay (sec)')
        plt.ylabel('spatial shift (PIV box unit)')
        plt.plot(tab_times_after,(opt[0]+np.sqrt(pcov[0,0]))*tab_times_after)
        plt.plot(tab_times_after,(opt[0]-np.sqrt(pcov[0,0]))*tab_times_after)
        plt.show()
    return opt,pcov


def convert_velocity_from_pxsec_to_mps(v_px_sec,dcm,dpx,w_piv=16):
    return v_px_sec * (w_piv/2) * (dcm*1e-2)/dpx#factor_px_to_meters