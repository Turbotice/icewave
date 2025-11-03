# -*- coding: utf-8 -*-
"""
Created on Mon Jun 23 15:08:13 2025

@author: sebas
"""

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
import scipy

import icewave.tools.matlab_colormaps as matcmaps
parula_map = matcmaps.parula()

global g 
g = 9.81


def dim_numbers(params,g = 9.81):
    """ Return dimensionless numbers from set of parameters
    Inputs: - params, tuple, (f_ex,h,rho_1,rho_2,nu_1,nu_2)
    Output: - dimensionless, tuple, (gamma,delta_1,delta_2,r) """
    
    omega = 2*np.pi*params[0]
    k_0 = omega**2/g
    
    rho_1 = params[2]
    rho_2 = params[3]
    nu_1 = params[4]
    nu_2 = params[5]
    
    gamma = params[1]*k_0
    delta_1 = np.sqrt(nu_1/omega)*k_0
    delta_2 = np.sqrt(nu_2/omega)*k_0
    r = rho_1/rho_2
    dimensionless = (gamma,delta_1,delta_2,r)
    return dimensionless

#----------------------------------------------------------------------------------------------------------------

def define_M(kappa,dimensionless):
    """ Create matrix for dispersion relation 
    Inputs : - k, float, wavevector value
             - dimensionless, parameters used to build matrix M (gamma,delta_1,delta_2,r)
             
    Outputs : - M, numpy.array, matrix of coefficients used to get dispersion relation. 
        This matrix construction is based on a vector [A_1,B_1,C_1,D_1,A_2,C_2] """

    
    # recover dimensionless parameters
    gamma = dimensionless[0]
    delta_1 = dimensionless[1]
    delta_2 = dimensionless[2]
    r = dimensionless[3]
    
    # define parameters for matrix M
    m_1 = np.sqrt(kappa**2 - 1j*1/delta_1**2)
    m_2 = np.sqrt(kappa**2 - 1j*1/delta_2**2)
    
    alpha_1 = delta_1*kappa
    alpha_2 = delta_2*kappa
    
    K_1 = m_1/kappa
    K_2 = m_2/kappa
    
    # define Matrix M
    # typical vector is (A1,B1,C1,D1,A2,C2)
      
    M = np.array([
        [(2*alpha_1**2*K_1 + 1j*kappa), 
         (-2*alpha_1**2*K_1 + 1j*kappa), 
         -1 + kappa - 1j*2*alpha_1**2, 
         -1 - kappa - 2*1j*alpha_1**2,
         0,
         0], # row 1
        
        [1j*alpha_1**2*(1 + K_1**2), 
         1j*alpha_1**2*(1 + K_1**2), 
         2*alpha_1**2, 
         -2*alpha_1**2,
         0,
         0], # row 2
        
        [1j*K_1*np.exp(-m_1*gamma), 
         -1j*K_1*np.exp(m_1*gamma), 
         np.exp(-kappa*gamma), 
         np.exp(+kappa*gamma),
         -1j*r*K_2*np.exp(-m_2*gamma),
         -r*np.exp(-kappa*gamma)], #row 3
        
        [1j*np.exp(-m_1*gamma), 
         1j*np.exp(m_1*gamma),
         np.exp(-kappa*gamma),
         -np.exp(+kappa*gamma),
         -1j*r*np.exp(-m_2*gamma),
         -r*np.exp(-kappa*gamma)], # row 4 
        
        [-2*alpha_1**2*K_1*np.exp(-m_1*gamma), 
         +2*alpha_1**2*K_1*np.exp(m_1*gamma),
         (1 + 1j*2*alpha_1**2)*np.exp(-kappa*gamma), 
         (1 + 1j*2*alpha_1**2)*np.exp(kappa*gamma),
         (2*alpha_2**2*K_2 + 1j*(1-r)*kappa)*np.exp(-m_2*gamma),
         (-1 - 1j*2*alpha_2**2 + (1-r)*kappa)*np.exp(-kappa*gamma)], # row 5
        
        [1j*alpha_1**2*(1 + K_1**2)*np.exp(-m_1*gamma), 
         1j*alpha_1**2*(1 + K_1**2)*np.exp(+m_1*gamma),
         2*alpha_1**2*np.exp(-kappa*gamma),
         -2*alpha_1**2*np.exp(+kappa*gamma),
         -1j*alpha_2**2*(1+K_2**2)*np.exp(-m_2*gamma),
         -2*alpha_2**2*np.exp(-kappa*gamma)] # row 6
        
        ])

    return M 

#-------------------------------------------------------------------------------------------------------------

def det_M(k,params):
    
    M = define_M(k,params)
    determinant = np.linalg.det(M)
    return determinant

#---------------------------------------------------------------------------------------------------------------

def d_det_M(k, params, h=1e-6*(1+1j)):
    # Numerical derivative using central difference
    return (det_M(k + h, params) - det_M(k - h, params)) / (2*h)

#---------------------------------------------------------------------------------------------------------------

def newton_raphson(params, k0, tol=1e-3, max_iter=100, step_derivative = 1e-6,threshold_derivative = 1e-12):
    
    """ Perform Newton-Raphson method. Enable to find value of k (complex) for which
    det(M(k)) = 0. 
    Inputs : - params, tuple, set of parameters used to define matrix M
             - k0, complex value, initial guess of k 
             - tol, float, tolerence below which we stop Newton algorithm
             - max_iter, integer, maximum number of iterations
             - step_derivative, float, step of wavevector for which we compute derivative
             - threshold_derivative, float, threshold below which the derivative is considered to be too small """
    
    k = k0
    for i in range(max_iter):
        f = det_M(k, params)
        df = d_det_M(k, params, h = step_derivative)
        # if abs(df) < threshold_derivative:
        #     raise ValueError("Derivative too small; might not converge.")
        k_new = k - f / df
        if abs(k_new - k) < tol:
            return k_new
        k = k_new
    raise RuntimeError("Newton-Raphson did not converge")
    

#---------------------------------------------------------------------------------------------------------------
def plot_det_M(x,y,dimensionless,cmap = parula_map):
    """ Compute determinant of matrix M for different values of k and plot the determinant in the complex plane
    Inputs : - x, numpy array, real part of wavevector
             - y, numpy array, imaginary part of wavevector
             - dimensionless, tuple (gamma,delta_1,delta_2,r)
             - cmap, colormap, optional colormap used to plot determinant
    Outputs : - fig, matplotlib figure
              - ax, matplotlib axis """
        
    X, Y = np.meshgrid(x, y)
    Z = X + 1j * Y

    # Evaluate |det(M(k))| over the grid

    det_values = np.zeros_like(Z, dtype=complex)

    for i in range(Z.shape[0]):
        for j in range(Z.shape[1]):
            k = Z[i, j]

            try:
                det_values[i, j] = det_M(k,dimensionless)
            except np.linalg.LinAlgError:
                det_values[i, j] = np.nan  # handle singularities or bad evals
                
    fig, ax = plt.subplots()
    extents = np.array((x.min(),x.max(),y.min(),y.max()))
    c = ax.imshow(np.log10(abs(det_values)),cmap = cmap,origin = 'lower',aspect = 'equal',
                  interpolation = 'gaussian',extent = extents)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2%", pad=0.1)
    
    cbar = plt.colorbar(c,cax = cax)
    cbar.set_label(r'$\mathrm{log10(|det(M(k))|)}$',labelpad = 5)
    ax.set_xlabel(r'$\mathrm{Re}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    ax.set_ylabel(r'$\mathrm{Im}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
    
    return fig, ax

#-------------------------------------------------------------------------------------------------------

def func_real_imag(z,dimensionless):
    """ Compute determinant of matrix M, for a given complex wavevector z and a set of dimensionless numbers 
    Inputs : - z, complex, wavector, complex
             - dimensionless, tuple (gamma,delta_1,delta_2,r)"""
    
    x,y = z # z = [Re(z),Im(z)]
    
    kappa = x + 1j*y
    determinant = det_M(kappa,dimensionless)
    return [np.real(determinant),np.imag(determinant)] # system of 2 real equations 

#-------------------------------------------------------------------------------------------------------

def get_wavevector_theory(freq,h,rho_1,rho_2,nu_1,nu_2):
    """ Compute wavector using viscous bilayers, for a given set of frequencies and 
    fluid physical properties. """
    
    params_array = np.array([(f_ex,h,rho_1,rho_2,nu_1,nu_2) for f_ex in freq])
    dimensionless_array = np.array([dim_numbers(params) for params in params_array]) 
    omega_array = 2*np.pi*freq
    k0_array = omega_array**2/g
    
    initial_guess = [1 , 0.1] # initial guess for [Re(kappa), Im(kappa)]
    k_array = np.zeros(len(freq),dtype = 'complex')
    for i,dimensionless in enumerate(dimensionless_array): # loop over dimensionless parameters
        
        solution = scipy.optimize.fsolve(func_real_imag,initial_guess,(dimensionless,))
        root = solution[0] + 1j*solution[1]
        print(f'Root of det(M) = {root}')
        initial_guess = solution
        k_array[i] = root*k0_array[i]
        
    return k_array


#--------------------------------------------------------------------------------------------------------

def define_M_hat(kappa,dimensionless):
    """ Create matrix for dispersion relation 
    Inputs : - kappa, float, dimensionless wavevector value
             - dimensionless, tuple of parameters used to build matrix M (gamma,delta_1,delta_2,r)
             
    Outputs : - M, numpy.array, matrix of coefficients used to get dispersion relation. 
        This matrix construction is based on a vector [A_1,B_1,C_1,D_1,A_2,C_2] """

   
    # recover dimensionless parameters
    gamma = dimensionless[0]
    delta_1 = dimensionless[1]
    delta_2 = dimensionless[2]
    r = dimensionless[3]
    
    # define parameters for matrix M
    m_1 = np.sqrt(kappa**2 - 1j*1/delta_1**2)
    m_2 = np.sqrt(kappa**2 - 1j*1/delta_2**2)
    
    alpha_1 = delta_1*kappa
    alpha_2 = delta_2*kappa
    
    K_1 = m_1/kappa
    K_2 = m_2/kappa
    # define Matrix M
    # typical vector is (A1,B1,C1,D1,A2,C2)
      
    M = np.array([
        [2*alpha_1**2*K_1 , 
         -2*alpha_1**2*K_1, 
         -1 - 1j*2*alpha_1**2, 
         -1 - 2*1j*alpha_1**2,
         0,
         0], # row 1
         
         
         # [1j*rho_1*nu_1, 
         #  1j*rho_1*nu_1,
         #  +1,
         #  -1,
         #  0,
         #  0], # row 1
        
        [1j*alpha_1**2*(1 + K_1**2),
         1j*alpha_1**2*(1 + K_1**2),
         2*alpha_1**2,
         -2*alpha_1**2,
         0,
         0], # row 2
        
        [1j*K_1*np.exp(-m_1*gamma),
         -1j*K_1*np.exp(m_1*gamma), 
         np.exp(-kappa*gamma), 
         np.exp(+kappa*gamma),
         -1j*r*K_2*np.exp(-m_2*gamma),
         -r*np.exp(-kappa*gamma)], #row 3
        
        [1j*np.exp(-m_1*gamma), 
         1j*np.exp(m_1*gamma), 
         +np.exp(-kappa*gamma), 
         -np.exp(+kappa*gamma),
         -1j*r*np.exp(-m_2*gamma),
         -r*np.exp(-kappa*gamma)], # row 4 
        
        [-2*alpha_1**2*K_1*np.exp(-m_1*gamma), 
         2*alpha_1**2*K_1*np.exp(m_1*gamma),
         (1 + 1j*2*alpha_1**2)*np.exp(-kappa*gamma),
         (1 + 1j*2*alpha_1**2)*np.exp(kappa*gamma),
         (2*alpha_2**2*K_2 + 1j*(1-r)*kappa)*np.exp(-m_2*gamma),
         (-1 - 1j*2*alpha_2**2 + (1-r)*kappa)*np.exp(-kappa*gamma)], # row 5
        
        [1j*alpha_1**2*(1+K_1**2)*np.exp(-m_1*gamma),
         1j*alpha_1**2*(1 + K_1**2)*np.exp(+m_1*gamma),
         2*alpha_1**2*np.exp(-kappa*gamma),
         -2*alpha_1**2*np.exp(+kappa*gamma),
         -1j*alpha_2**2*(1 + K_2**2)*np.exp(-m_2*gamma),
         -2*alpha_2**2*np.exp(-kappa*gamma)] # row 6
        
        ])

    return M 

#----------------------------------------------------------------------------------------------------------------

def get_coeffs(kappa,dimensionless):
    """ Get coefficients A1,B1,C1,D1,A2,C2 used to define pressure and vorticity fields.
    This methods solver a system MX = Y, based on the knowledge of the surface amplitude xi_0 and assuming that 
    surface and fluids interface deformation are in phase. 
    Inputs : - xi_0, float, surface amplitude in meters
             - k, complex, wavevector in rad/m
             - dimensionless, tuple, (gamma,delta_1,delta_2,r) """
             
    # Compute matrix for determination of coefficients
    M_hat = define_M_hat(kappa, dimensionless)
    
    Y = np.array([-1,0,0,0,0,0])
    # Y = np.array([rho_1*omega**2/new_k*xi_0,0,0,0,0,0]) 
    
    # Solve the system M * X = Y
    try:
        coeffs = np.linalg.solve(M_hat, Y)
        print("Solution X:", coeffs)
    except np.linalg.LinAlgError as e:
        print("Error:", e)
        
    return coeffs 
