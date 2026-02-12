# -*- coding: utf-8 -*-
"""
Created on Fri Jun 20 12:00:47 2025

@author: sebas
"""

import numpy as np

global g
g = 9.81

def define_M(k,params):
    """ Create matrix for dispersion relation 
    Inputs : - k, float, wavevector value
             - params, parameters used to build matrix M [f_ex,h,rho_1,rho_2,nu_1,nu_2]
             
    Outputs : - M, numpy.array, matrix of coefficients used to get dispersion relation. 
        This matrix construction is based on a vector [A_1,B_1,C_1,D_1,A_2,C_2] """

    f_ex = params[0]
    h = params[1]
    rho_1 = params[2]
    rho_2 = params[3]
    nu_1 = params[4]
    nu_2 = params[5]             
            
    omega = 2*np.pi*f_ex # in rad/s
    
    sigma = g*k/omega**2
    r = rho_1/rho_2
    
    # define dimensionless parameters
    m_1 = np.sqrt(k**2 - 1j*omega/nu_1)
    m_2 = np.sqrt(k**2 - 1j*omega/nu_2)
    
    delta_1 = np.sqrt(nu_1/omega)
    delta_2 = np.sqrt(nu_2/omega)
    
    K_1 = m_1/k
    K_2 = m_2/k
    
    alpha_1 = delta_1*k
    alpha_2 = delta_2*k
    
    gamma = k*h
    
    # define Matrix M
    # typical vector is (A1,B1,C1,D1,A2,C2)
      
    M = np.array([
        [2*alpha_1**2*K_1 + 1j*sigma, 
         -2*alpha_1**2*K_1 + 1j*sigma, 
         -1 + sigma - 1j*2*alpha_1**2, 
         -1 - sigma - 2*1j*alpha_1**2,
         0,
         0], # row 1
        
        [1j*alpha_1**2*(1 + K_1**2), 
         1j*alpha_1**2*(1 + K_1**2), 
         2*alpha_1**2, 
         -2*alpha_1**2,
         0,
         0], # row 2
        
        [1j*K_1*np.exp(-K_1*gamma), 
         -1j*K_1*np.exp(K_1*gamma), 
         np.exp(-gamma), 
         np.exp(+gamma),
         -1j*r*K_2*np.exp(-K_2*gamma),
         -r*np.exp(-gamma)], #row 3
        
        [1j*np.exp(-K_1*gamma), 
         1j*np.exp(K_1*gamma),
         np.exp(-gamma),
         -np.exp(+gamma),
         -1j*r*np.exp(-K_2*gamma),
         -r*np.exp(-gamma)], # row 4 
        
        [-2*alpha_1**2*K_1*np.exp(-K_1*gamma), 
         +2*alpha_1**2*K_1*np.exp(K_1*gamma),
         (1 + 1j*2*alpha_1**2)*np.exp(-gamma), 
         (1 + 1j*2*alpha_1**2)*np.exp(gamma),
         (2*alpha_2**2*K_2 + 1j*(1-r)*sigma)*np.exp(-K_2*gamma),
         (-1 - 1j*2*alpha_2**2 + (1-r)*sigma)*np.exp(-gamma)], # row 5
        
        [1j*rho_1*nu_1*alpha_1**2*(1 + K_1**2)*np.exp(-K_1*gamma), 
         1j*rho_1*nu_1*alpha_1**2*(1 + K_1**2)*np.exp(+K_1*gamma),
         2*alpha_1**2*np.exp(-gamma),
         -2*alpha_1**2*np.exp(+gamma),
         -1j*alpha_2**2*(1+K_2**2)*np.exp(-K_2*gamma),
         -2*alpha_2**2*np.exp(-gamma)] # row 6
        
        ])

    return M 


def define_M_hat(k,params):
    
    f_ex = params[0]
    h = params[1]
    rho_1 = params[2]
    rho_2 = params[3]
    nu_1 = params[4]
    nu_2 = params[5]             
            
    omega = 2*np.pi*f_ex # in rad/s
    
    sigma = g*k/omega**2
    r = rho_1/rho_2
    
    # define dimensionless parameters
    m_1 = np.sqrt(k**2 - 1j*omega/nu_1)
    m_2 = np.sqrt(k**2 - 1j*omega/nu_2)
    
    delta_1 = np.sqrt(nu_1/omega)
    delta_2 = np.sqrt(nu_2/omega)
    
    K_1 = m_1/k
    K_2 = m_2/k
    
    alpha_1 = delta_1*k
    alpha_2 = delta_2*k
    
    gamma = k*h
    # define Matrix M
    # typical vector is (A1,B1,C1,D1,A2,C2)
      
    M = np.array([
        [2*rho_1*nu_1*alpha_1**2*K_1 , 
         -2*rho_1*nu_1*alpha_1**2*K_1, 
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
        
        [1j*rho_1*nu_1*alpha_1**2*(1 + K_1**2),
         1j*rho_1*nu_1*alpha_1**2*(1 + K_1**2),
         2*alpha_1**2,
         -2*alpha_1**2,
         0,
         0], # row 2
        
        [1j*rho_1*nu_1*K_1*np.exp(-K_1*gamma),
         -1j*rho_1*nu_1*K_1*np.exp(K_1*gamma), 
         np.exp(-gamma), 
         np.exp(+gamma),
         -1j*r*rho_2*nu_2*K_2*np.exp(-K_2*gamma),
         -r*np.exp(-gamma)], #row 3
        
        [1j*rho_1*nu_1*np.exp(-K_1*gamma), 
         1j*rho_1*nu_1*np.exp(K_1*gamma), 
         +np.exp(-gamma), 
         -np.exp(+gamma),
         -1j*r*rho_2*nu_2*np.exp(-K_2*gamma),
         -r*np.exp(-gamma)], # row 4 
        
        [-2*rho_1*nu_1*alpha_1**2*K_1*np.exp(-K_1*gamma), 
         2*rho_1*nu_1*alpha_1**2*K_1*np.exp(K_1*gamma),
         (1 + 1j*2*alpha_1**2)*np.exp(-gamma),
         (1 + 1j*2*alpha_1**2)*np.exp(gamma),
         rho_2*nu_2*(2*alpha_2**2*K_2 + 1j*(1-r)*sigma)*np.exp(-K_2*gamma),
         (-1 - 1j*2*alpha_2**2 + (1-r)*sigma)*np.exp(-gamma)], # row 5
        
        [1j*rho_1*nu_1*alpha_1**2*(1+K_1**2)*np.exp(-K_1*gamma),
         1j*rho_1*nu_1*alpha_1**2*(1 + K_1**2)*np.exp(+K_1*gamma),
         2*alpha_1**2*np.exp(-gamma),
         -2*alpha_1**2*np.exp(+gamma),
         -1j*rho_2*nu_2*alpha_2**2*(1 + K_2**2)*np.exp(-K_2*gamma),
         -2*alpha_2**2*np.exp(-gamma)] # row 6
        
        ])

    return M 

def det_M(k,params):
    
    M = define_M(k,params)
    determinant = np.linalg.det(M)
    return determinant

def d_det_M(k, params, h=1e-6*(1+1j)):
    # Numerical derivative using central difference
    return (det_M(k + h, params) - det_M(k - h, params)) / (2*h)

def newton_raphson(params, k0, tol=1e-8, max_iter=100, step_derivative = 1e-6,threshold_derivative = 1e-12):
    
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
        if abs(df) < threshold_derivative:
            raise ValueError("Derivative too small; might not converge.")
        k_new = k - f / df
        if abs(k_new - k) < tol:
            return k_new
        k = k_new
    raise RuntimeError("Newton-Raphson did not converge")
    
def get_coeffs(xi_0,k,params):
    """ Get coefficients A1,B1,C1,D1,A2,C2 used to define pressure and vorticity fields.
    This methods solver a system MX = Y, based on the knowledge of the surface amplitude xi_0 and assuming that 
    surface and fluids interface deformation are in phase. 
    Inputs : - xi_0, float, surface amplitude in meters
             - k, complex, wavevector in rad/m
             - params, tuple, (f_ex,h,rho_1,rho_2,nu_1,nu_2) """
             
    # Compute matrix for determination of coefficients
    M_hat = define_M_hat(k, params)
    
    Y = np.array([-params[2]*g*xi_0,0,0,0,0,0])
    # Y = np.array([rho_1*omega**2/new_k*xi_0,0,0,0,0,0]) 
    
    # Solve the system M * X = Y
    try:
        coeffs = np.linalg.solve(M_hat, Y)
        print("Solution X:", coeffs)
    except np.linalg.LinAlgError as e:
        print("Error:", e)
        
    return coeffs  