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

# %% fonctions pour importer les données matlab
def extract_data_mat_file(dd):
    with h5py.File(dd['path_mat_file'], 'r') as file:
        data = file['data'][:]  # Accessing 'data' variable
    dd['data'] = data
#    return data

def tuples_to_complex(dd):
    matrix = dd['data']
    # Initialize an empty matrix to store complex numbers
    complex_matrix = np.empty(matrix.shape, dtype=complex)
    
    # Iterate through each element in the matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            # Get the tuple from the original matrix
            a, b = matrix[i, j]
            
            # Create a complex number from the tuple and assign it to the corresponding position
            complex_matrix[i, j] = complex(a, b)
    dd['complex_field'] = complex_matrix    
    dd['real_field_t0'] = np.real(complex_matrix)/np.abs(complex_matrix) # peut être utile

def compute_real_field(dd,t):
    complex_field = dd['complex_field']
    freq = dd['f_exc']
    return np.real(complex_field * np.exp(-1j*2*np.pi*freq*t))/np.abs(complex_field)

#    return complex_matrix

# %% fonction pour affecter tous les parametres dans des variables globales plus pratiques
def affect_params(dd):
    global smooth
    global window_size
    global order
    global complex_field
    global freq
    global coord_source
    global f_exc
    global T
    global unwrap
    global shift_startpoint
    global plots
    global Vmax
    global index_point_on_circle
    smooth = dd['savgol_param']['smooth']
    window_size = dd['savgol_param']['window_size']
    order = dd['savgol_param']['order']
    complex_field = dd['complex_field']
    freq = dd ['f_exc']
    f_exc = dd['f_exc']
    T = 1/freq
    dd['T'] = T
    unwrap = dd['unwrap']
    shift_startpoint = dd['shift_startpoint']
    plots = dd['plots']
    Vmax = dd['Vmax']
    index_point_on_circle = dd['index_point_on_circle']
    #print(Vmax)
    if ('x_circle_on_image' in dd)&('y_circle_on_image' in dd):
        global x_circle_on_image
        global y_circle_on_image
        x_circle_on_image = dd['x_circle_on_image']
        y_circle_on_image = dd['y_circle_on_image']
        coord_source = dd['coord_source']

# %% fonctions pour tracer cercle sur image etc.
def Circle(dd):#coord de O dans l'ordre x_coord,y_coord (où l'axe y est inversé)
    xvals = np.arange(dd['complex_field'].shape[1])
    O = dd['coord_source']
    R = dd['radius']
    yp = O[1] + np.sqrt(R**2 - (xvals-O[0])**2)
    ym = O[1] - np.sqrt(R**2 - (xvals-O[0])**2)
    x_circle = np.hstack((xvals,xvals))
    y_circle = np.hstack((yp,ym))
    return x_circle,y_circle

def coord_circle_on_image(dd):
    affect_params(dd)
    real_field_t0 = np.real(complex_field)/np.abs(complex_field)
    x_circle,y_circle = Circle(dd)
    indices_to_plot = (x_circle>0)&(x_circle<complex_field.shape[1])&(y_circle>0)&(y_circle<complex_field.shape[0])

    x_circle_on_image = x_circle[indices_to_plot]
    y_circle_on_image = y_circle[indices_to_plot]

    if plots==True:
        plt.figure()
        plt.imshow(real_field_t0)
        plt.plot(x_circle_on_image,y_circle_on_image,'r')
        plt.xlabel('x',fontsize=15)
        plt.ylabel('y',fontsize=15)
        #plt.ylim(-50,120)
        plt.show()
    
    dd['x_circle_on_image'] = x_circle_on_image
    dd['y_circle_on_image'] = y_circle_on_image

# %% fonctions pour calculer les profils du champ de hauteur du point source à un point donné sur le cercle##
def field_along_profile(dd,coord_start_point,coord_end_point,t):
    image = compute_real_field(dd,t)
    if smooth==True:
        window_size_new = int(window_size/2) * 2 + 1
        image = sgolay2d(image,window_size_new,order)
    image_with_nans = np.vstack((np.ones_like(image)*np.nan , image))

    coord_start_point_with_nans = (coord_start_point[0],coord_start_point[1]+image.shape[0])
    coord_end_point_with_nans = (coord_end_point[0],coord_end_point[1]+image.shape[0])

    p = profile_line(image_with_nans,np.flip(coord_start_point_with_nans),np.flip(coord_end_point_with_nans),cval=np.nan)
    
    if plots==True:
        plt.figure()
        plt.imshow(image_with_nans,vmin=-1,vmax=1)
        plt.plot(coord_start_point_with_nans[0],coord_start_point_with_nans[1],'ob')
        plt.plot(coord_end_point_with_nans[0],coord_end_point_with_nans[1],'or')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.colorbar()
        plt.show()        
    return p

def field_along_profile_endpointoncircle(dd,index_point_on_circle,t): # attention le probleme est que ca ne doit pas depasser la longueur du tableau x_circle_on_image
    affect_params(dd)   
    coord_start_point = coord_source
    coord_end_point = (np.round(x_circle_on_image[index_point_on_circle]).astype(int),np.round(y_circle_on_image[index_point_on_circle]).astype(int))
    p = field_along_profile(dd,coord_start_point,coord_end_point,t)
    return p

# %% savgol smoothing function
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


# %% fonctions pour correler les profils et en déduire un décalage spatial entre les deux profils (marche moyennement bien)
    
def correlate_2profiles(dd,time_before,time_after,index_point_on_circle):
    affect_params(dd)
    image_a = compute_real_field(dd,time_before)
    image_b = compute_real_field(dd,time_after)
#    if smooth==True:
#        window_size_new = int(window_size/2) * 2 + 1
#        image_a = sgolay2d(image_a,window_size_new,order)
#        image_b = sgolay2d(image_b,window_size_new,order)
    coord_source_new = coord_source
    x_circle_on_image_new = x_circle_on_image
    y_circle_on_image_new = y_circle_on_image
    p_a = field_along_profile_endpointoncircle(dd,index_point_on_circle,time_before)[shift_startpoint:]
    p_b = field_along_profile_endpointoncircle(dd,index_point_on_circle,time_after)[shift_startpoint:]
    p_a_with_zeros = np.where(np.isnan(p_a),0,p_a)
    p_b_with_zeros = np.where(np.isnan(p_b),0,p_b)
    p_a_on_image = np.delete(p_a_with_zeros,np.where(np.isnan(p_a))) ####### ATTENTION
    p_b_on_image = np.delete(p_b_with_zeros,np.where(np.isnan(p_b))) ####### ATTENTION
    correlation_a_b = signal.correlate(p_a_on_image,p_b_on_image,mode='full')
    lags = signal.correlation_lags(len(p_a_on_image), len(p_b_on_image)) # voir comment sont calculés les lags
    return lags,correlation_a_b,p_a_on_image,p_b_on_image

def compute_dist_max_corr(dd,time_before,time_after,index_point_on_circle):
    affect_params(dd)

    lags,correlation_a_b,p_a_on_image,p_b_on_image = correlate_2profiles(dd,time_before,time_after,index_point_on_circle)
    index_max_correlation_a_b = np.where(correlation_a_b==np.max(correlation_a_b))
    dist_max_corr = lags[index_max_correlation_a_b] # distance par unité de cases PIV
    if plots==True:
        plt.figure()
        plt.plot(p_a_on_image)
        plt.plot(p_b_on_image)
        #plt.show()
        plt.figure()
        plt.plot(lags,correlation_a_b)
        print(dist_max_corr)
        plt.vlines(dist_max_corr,np.min(correlation_a_b),np.max(correlation_a_b),linewidth=0.7,linestyle='--',color='red')
        plt.show()
    return dist_max_corr,correlation_a_b,p_a_on_image,p_b_on_image

# fit qui utilise plusieurs dist_max_corr par rapport à l'instant 0
def linear(x,a):
    return a*x

def correlate_full_period_with_t0(dd,tab_times_after):
    affect_params(dd)
    tab_dist_max_corr = np.zeros_like(tab_times_after)
    for i in range(len(tab_times_after)):
        dist_max_corr = compute_dist_max_corr(complex_field,0,tab_times_after[i],index_point_on_circle,coord_source,x_circle_on_image,y_circle_on_image,freq,shift_startpoint=shift_startpoint,plot=plot,smooth=smooth)[0]
        #print(dist_max_corr[0])
        tab_dist_max_corr[i] = dist_max_corr[0]
    return tab_dist_max_corr

# %% define a function that detects the wave front with subpixelar precision
def compute_lag_sub_pix(dd,a,b):
    #print(a)
    #print(b)
    affect_params(dd)
    a = np.delete(a,np.where(np.isnan(a)))
    b = np.delete(b,np.where(np.isnan(b)))
    #print(Vmax)
    P1 = np.zeros(2*Vmax-1)
    P2 = np.zeros(2*Vmax-1)
    for v in range(-Vmax+1, Vmax):
#        print((a[Vmax-v-1:len(a)-Vmax-v+1] - b[Vmax-1:len(a)-Vmax+1])**2)
        P1[v+Vmax-1] = np.sum((a[Vmax-v-1:len(a)-Vmax-v+1] - b[Vmax-1:len(a)-Vmax+1])**2)
        P2[v+Vmax-1] = np.sum((b[Vmax-v-1:len(a)-Vmax-v+1] - a[Vmax-1:len(a)-Vmax+1])**2)
    P2 = np.flip(P2)
    P = P1 + P2
    #print('length of P',len(P))
    #print('P:',P)

    i_min = np.argmin(P)
    delta_i = (P[i_min+1] - P[i_min-1]) / (2 * (2 * P[i_min] - P[i_min+1] - P[i_min-1]))
    lag = i_min - Vmax + delta_i + 1 # il faut rajouter 1 car argmin donne l'indice du minimum mais en convention python
    if plots==True:
        plt.figure()
        plt.plot(a[Vmax-v:len(a)-Vmax-v+1])
        plt.plot(b[Vmax-v:len(a)-Vmax-v+1])
        plt.show()
        plt.plot(np.arange(1, len(P)+1)-Vmax, P, '+')
        plt.show()


    return lag,P

def compute_profile_and_lag_sub_pix(dd,time_before,time_after):#,index_point_on_circle):
    affect_params(dd)
    p_a = field_along_profile_endpointoncircle(dd,index_point_on_circle,time_before)
    p_b = field_along_profile_endpointoncircle(dd,index_point_on_circle,time_after)
    lag,P = compute_lag_sub_pix(dd,p_a,p_b)
    return lag,P

# fit qui utilise plusieurs dist_max_corr par rapport à l'instant 0

def compute_lags_subpix_time_interval_with_t0(dd,t0,tab_times_after):#,index_point_on_circle):
    affect_params(dd)
    tab_lags = np.zeros_like(tab_times_after)
    for i in range(len(tab_times_after)):
        tab_lags[i] = compute_profile_and_lag_sub_pix(dd,t0,tab_times_after[i])[0]
    return tab_lags

def compute_LAGS_subpix_full_period(dd,nb_intervals):
    affect_params(dd)
    width_time_interval = T/nb_intervals
    LAGS = []
    TAB_TIMES_AFTER = []
    for i in range(nb_intervals):
        t0 = i*width_time_interval
        tab_times_after = np.linspace(t0,t0+width_time_interval,5)
        tab_lags = compute_lags_subpix_time_interval_with_t0(dd,t0,tab_times_after)
        LAGS.append(tab_lags)
        TAB_TIMES_AFTER.append(tab_times_after)
    return LAGS,TAB_TIMES_AFTER
# %% use the above function to compute the phase velocity (for the correlation method)
def compute_phase_velocity(tab_times_after,tab_dist_max_corr):
    affect_params(dd)
    if unwrap==True:
        tab_dist_max_corr = np.unwrap(tab_dist_max_corr)
    #dt = tab_times_after[1] - tab_times_after[0]
    opt,pcov = curve_fit(linear,tab_times_after,tab_dist_max_corr)
    if plots==True:
        plt.figure()
        plt.plot(tab_times_after,tab_dist_max_corr,'+')
        plt.plot(tab_times_after,linear(tab_times_after,opt[0]))
        plt.xlabel('delay (sec)')
        plt.ylabel('spatial shift (PIV box unit)')
        plt.plot(tab_times_after,(opt[0]+np.sqrt(pcov[0,0]))*tab_times_after)
        plt.plot(tab_times_after,(opt[0]-np.sqrt(pcov[0,0]))*tab_times_after)
        plt.show()
    return opt,pcov


def convert_velocity_from_pivboxsec_to_mps(dd,v_pivbox_sec):
    w_piv = dd['W']
    dcm = dd['dcm']
    dpx = dd['dpx']
    return v_pivbox_sec * (w_piv/2) * (dcm*1e-2)/dpx#factor_px_to_meters

# %% use the above function to compute the phase velocity (for the subpix method)
def compute_phase_velocity_sub_pix(dd):
    affect_params(dd)
    LAGS = compute_LAGS_subpix_full_period(dd,nb_intervals=10)[0]
    #for i in range(len(LAGS)):
    #    width_time_interval = T/len(LAGS[i])
    #    t0 = i*width_time_interval
    #    tab_times_after = np.linspace(t0,t0+width_time_interval,nb_points_per_interval)
    #    tab_lags = LAGS[i]
        #plt.plot(tab_times_after,tab_lags,'^')
    t_fit = np.repeat(np.linspace(0,width_time_interval,nb_points_per_interval),nb_time_intervals)
    lags_fit = (np.array(LAGS).T).flatten()

    opt,pcov = curve_fit(linear,t_fit,lags_fit)
    v_phase = convert_velocity_from_pivboxsec_to_mps(dd,opt[0])
    u_v_phase = convert_velocity_from_pivboxsec_to_mps(dd,np.sqrt(pcov[0][0]))
    if plots==True:
        plt.figure()
        plt.plot(t_fit,linear(t_fit,opt),'r')
        plt.plot(t_fit,lags_fit,'o')
        plt.plot(t_fit,linear(t_fit,opt+np.sqrt(pcov)[0]),'k--')
        plt.plot(t_fit,linear(t_fit,opt-np.sqrt(pcov)[0]),'k--')
        plt.show()

    return LAGS, v_phase, u_v_phase, t_fit, lags_fit

# %% def params
dd = {}
dd['f_exc'] = 100#150#
#dd['T'] = 1/dd['f_exc'] # ne pas executer cette ligne car T est calculé directement lorqu'on affecte f_exc normalement
dd['freq_acq'] = 99.0099#148.5222#
dd['dcm'] = 7
dd['dpx'] = 1192
dd['W'] = 32
dd['Dt'] = 8
#dd['date'] = '20240419'
dd['date'] = '20240514'
dd['main_folder'] = 'X:/Banquise/Vasco/Frigo_pmmh/'+dd['date']+'/'
dd['path_mat_file'] = dd['main_folder']+str(dd['f_exc'])+'Hz_'+str(dd['freq_acq'])+'Hz/matData/video_demod_W'+str(dd['W'])+'_Dt'+str(dd['Dt'])+'/figdata_complex.mat'
dd['coord_source'] = (94,-111)#(70,-130)
dd['radius'] = 235
dd['savgol_param'] = {'smooth':True,'window_size':41,'order':4}
dd['plots'] = True
dd['shift_startpoint'] = 0
dd['unwrap'] = True
dd['Vmax'] = 20
dd['index_point_on_circle'] = 20
# %% read data

extract_data_mat_file(dd)
tuples_to_complex(dd)
Circle(dd)
coord_circle_on_image(dd)
affect_params(dd)

# %% test
#affect_params(dd)
#dist_max_corr,correlation_a_b,p_a_on_image,p_b_on_image = compute_dist_max_corr(dd,0,0.1*T,30)
lag,P = compute_profile_and_lag_sub_pix(dd,0,0.1*T)


# %%
tab_f_exc = np.array([50,70,100,120,150,170])
tab_freq_acq = np.array([49.5049,69.30487,99.0099,118.8072,148.5222,84.5738])

tab_v_phase = np.zeros(len(tab_f_exc))

for i in range(len(tab_f_exc)):
    dd['f_exc'] = tab_f_exc[i]
    dd['freq_acq'] = tab_freq_acq[i]
    dd['plots']=False
    dd['index_point_on_circle'] = 20
    LAGS, v_phase, u_v_phase, t_fit, lags_fit = compute_phase_velocity_sub_pix(dd)
    #v_phase = v_phase/2 # A ENLEVER UNE FOIS QUE LES DONNES SERONT CORRIGEES POUR LA VALEUR DE W !!!!!
    tab_v_phase[i] = v_phase


 #%%
plt.figure()
plt.plot(tab_f_exc,tab_v_phase,'o')
#plt.loglog()
plt.show()
# %%
