# -*- coding: utf-8 -*-
"""
Created on Wed Jun 11 17:46:06 2025

@author: sebas
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import os

import icewave.tools.matlab_colormaps as matcmaps
import icewave.sebastien.set_graphs as set_graphs


# PARULA COLORMAP 
parula_map = matcmaps.parula()

plt.rc('text', usetex=True)
plt.rc('font', family='serif', serif='Computer Modern')

# Set global parameters
global g
g = 9.81

main_path = 'C:/Users/sebas/OneDrive/Bureau/These PMMH/Meetings/Meeting_Aurore_frasil/'
fig_folder = f'{main_path}model_2_layers_viscous_viscous/'
if not os.path.isdir(fig_folder):
    os.mkdir(fig_folder)
    
#%% Def function 

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
        [rho_1*nu_1*(2*alpha_1**2*K_1 + 1j*sigma), 
         rho_1*nu_1*(-2*alpha_1**2*K_1 + 1j*sigma), 
         -1 + sigma - 1j*2*alpha_1**2, 
         -1 - sigma - 2*1j*alpha_1**2,
         0,
         0], # row 1
        
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
         np.exp(-gamma),
         -np.exp(+gamma),
         -1j*r*rho_2*nu_2*np.exp(-K_2*gamma),
         -r*np.exp(-gamma)], # row 4 
        
        [-2*rho_1*nu_1*alpha_1**2*K_1*np.exp(-K_1*gamma), 
         +2*rho_1*nu_1*alpha_1**2*K_1*np.exp(K_1*gamma),
         (1 + 1j*2*alpha_1**2)*np.exp(-gamma), 
         (1 + 1j*2*alpha_1**2)*np.exp(gamma),
         rho_2*nu_2*(2*alpha_2**2*K_2 + 1j*(1-r)*sigma)*np.exp(-K_2*gamma),
         (-1 - 1j*2*alpha_2**2 + (1-r)*sigma)*np.exp(-gamma)], # row 5
        
        [1j*rho_1*nu_1*alpha_1**2*(1 + K_1**2)*np.exp(-K_1*gamma), 
         1j*rho_1*nu_1*alpha_1**2*(1 + K_1**2)*np.exp(+K_1*gamma),
         2*alpha_1**2*np.exp(-gamma),
         -2*alpha_1**2*np.exp(+gamma),
         -1j*rho_2*nu_2*alpha_2**2*(1+K_2**2)*np.exp(-K_2*gamma),
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

def compute_p1(z,k,coeffs):
    """ Compute pressure profile in fluid number 1
    Inputs : - z, array of vertical coordinates (zero at surface)
             - k, float, wavevector
             - coeffs, array of coefficients obtained from system resolution using wave amplitude as input, 
                 ordered as [A1,B1,C1,D1,A2,C2]
    Outputs : - p1(z), array of pressure vertical profile """
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    p1 = C1*np.exp(k*z) + D1*np.exp(-k*z)
    return p1
    
def compute_p2(z,k,coeffs):
    """ Compute pressure profile in fluid number 2
    Inputs : - z, array of vertical coordinates (zero at surface)
             - k, float, wavevector
             - coeffs, array of coefficients obtained from system resolution using wave amplitude as input
                 ordered as [A1,B1,C1,D1,A2,C2]
    Outputs : - p2(z), array of pressure vertical profile """    
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    p2 = C2*np.exp(k*z)
    return p2

def compute_vort1(z,params,k,coeffs):
    """ Compute vorticity in fluid number 1
    Inputs : - z, array of vertical coordinates (zero at surface)
             - params, array of physical parameters,[f_ex,h,rho_1,rho_2,nu_1,nu_2]
             - k, float, wavevector
             - coeffs, array of coefficients obtained from system resolution using wave amplitude as input
                 ordered as [A1,B1,C1,D1,A2,C2]
    Outputs : - vort1(z), array of vorticity vertical profile """

    [A1,B1,C1,D1,A2,C2] = coeffs
    f_ex = params[0]
    nu_1 = params[4]            
            
    omega = 2*np.pi*f_ex # in rad/s
    m_1 = np.sqrt(k**2 - 1j*omega/nu_1)
    
    vort1 = A1*np.exp(m_1*z) + B1*np.exp(-m_1*z)
    
    return vort1

def compute_vort2(z,params,k,coeffs):
    """ Compute vorticity in fluid number 2
    Inputs : - z, array of vertical coordinates (zero at surface)
             - params, array of physical parameters, [f_ex,h,rho_1,rho_2,nu_1,nu_2]
             - k, float, wavevector
             - coeffs, array of coefficients obtained from system resolution using wave amplitude as input
                 ordered as [A1,B1,C1,D1,A2,C2]
    Outputs : - vort2(z), array of vorticity vertical profile """
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    f_ex = params[0]
    nu_2 = params[5]             
            
    omega = 2*np.pi*f_ex # in rad/s
    m_2 = np.sqrt(k**2 - 1j*omega/nu_2)
    
    vort2 = A2*np.exp(m_2*z)
    return vort2

def compute_u1(z,params,k,coeffs):
    """ Compute horizontal velocity in fluid number 1
    Inputs : - z, array of vertical coordinates (zero at surface)
             - params, array of physical parameters, [f_ex,h,rho_1,rho_2,nu_1,nu_2]
             - k, float, wavevector
             - coeffs, array of coefficients obtained from system resolution using wave amplitude as input
                 ordered as [A1,B1,C1,D1,A2,C2]
    Outputs : - u1(z), array of horizontal velocity, vertical profile"""
    
    f_ex = params[0]
    rho_1 = params[2]
    nu_1 = params[4]          
            
    omega = 2*np.pi*f_ex # in rad/s
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    
    m_1 = np.sqrt(k**2 - 1j*omega/nu_1)
    
    pressure_term = (C1*np.exp(k*z) + D1*np.exp(-k*z))*k/omega/rho_1
    vorticity_term = (A1*np.exp(m_1*z) - B1*np.exp(-m_1*z))*1j*nu_1*m_1/omega
    
    u_1 = pressure_term + vorticity_term
    return u_1


def compute_w1(z,params,k,coeffs):
    """ Compute vertical velocity in fluid number 1
    Inputs : - z, array of vertical coordinates (zero at surface)
             - params, array of physical parameters, [f_ex,h,rho_1,rho_2,nu_1,nu_2]
             - k, float, wavevector
             - coeffs, array of coefficients obtained from system resolution using wave amplitude as input
                 ordered as [A1,B1,C1,D1,A2,C2]
    Outputs : - w1(z), array of vertical velocity, vertical profile"""
    
    f_ex = params[0]
    rho_1 = params[2]
    nu_1 = params[4]            
            
    omega = 2*np.pi*f_ex # in rad/s
    m_1 = np.sqrt(k**2 - 1j*omega/nu_1)
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    
    pressure_term = (C1*np.exp(k*z) - D1*np.exp(-k*z))*(-1j)*k/omega/rho_1
    vorticity_term = (A1*np.exp(m_1*z) + B1*np.exp(-m_1*z))*nu_1*k/omega
    w_1 = pressure_term + vorticity_term
    return w_1


def compute_u2(z,params,k,coeffs):
    """ Compute horizontal velocity in fluid number 2
    Inputs : - z, array of vertical coordinates (zero at surface)
             - params, array of physical parameters, [f_ex,h,rho_1,rho_2,nu_1,nu_2]
             - k, float, wavevector
             - coeffs, array of coefficients obtained from system resolution using wave amplitude as input
                 ordered as [A1,B1,C1,D1,A2,C2]
    Outputs : - u2(z), array of horizontal velocity, vertical profile"""
    
    f_ex = params[0]
    rho_2 = params[3]
    nu_2 = params[5]             
            
    omega = 2*np.pi*f_ex # in rad/s
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    
    m_2 = np.sqrt(k**2 - 1j*omega/nu_2)
    
    pressure_term = C2*np.exp(k*z)*k/omega/rho_2
    vorticity_term = A2*np.exp(m_2*z)*1j*nu_2*m_2/omega
    
    u_2 = pressure_term + vorticity_term
    return u_2


def compute_w2(z,params,k,coeffs):
    """ Compute vertical velocity in fluid number 2
    Inputs : - z, array of vertical coordinates (zero at surface)
             - params, array of physical parameters, [f_ex,h,rho_1,rho_2,nu_1,nu_2]
             - k, float, wavevector
             - coeffs, array of coefficients obtained from system resolution using wave amplitude as input
                 ordered as [A1,B1,C1,D1,A2,C2]
    Outputs : - w2(z), array of vertical velocity, vertical profile"""
    
    f_ex = params[0]
    rho_2 = params[3]
    nu_2 = params[5]             
            
    omega = 2*np.pi*f_ex # in rad/s
    m_2 = np.sqrt(k**2 - 1j*omega/nu_2)
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    
    pressure_term = C2*np.exp(k*z)*(-1j)*k/omega/rho_2
    vorticity_term = A2*np.exp(m_2*z)*nu_2*k/omega
    w_2 = pressure_term + vorticity_term
    
    return w_2

def compute_sigma_zz_1(z,params,k,coeffs):
    
    f_ex = params[0]
    rho_1 = params[2]
    nu_1 = params[4]            
            
    omega = 2*np.pi*f_ex # in rad/s
    m_1 = np.sqrt(k**2 - 1j*omega/nu_1)
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    
    sigma = g*k/omega**2
    delta_1 = np.sqrt(nu_1/omega)   
    K_1 = m_1/k    
    alpha_1 = delta_1*k
    
    sigma_zz = -rho_1*nu_1*(2*alpha_1**2 * K_1 + 1j*sigma)*np.exp(m_1*z)*A1
    + rho_1*nu_1*(2*alpha_1**2 * K_1 - 1j*sigma)*np.exp(-m_1*z)*B1
    + (1 - sigma + 2*1j*alpha_1**2)*np.exp(k*z)*C1
    + (1 + sigma + 2*1j*alpha_1**2)*np.exp(-k*z)*D1
        
    return sigma_zz
    
def compute_sigma_zz_2(z,params,k,coeffs):
    
    f_ex = params[0]
    rho_2 = params[3]
    nu_2 = params[5]            
            
    omega = 2*np.pi*f_ex # in rad/s
    m_2 = np.sqrt(k**2 - 1j*omega/nu_2)
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    
    sigma = g*k/omega**2
    delta_2 = np.sqrt(nu_2/omega)   
    K_2 = m_2/k    
    alpha_2 = delta_2*k
    
    sigma_zz = -rho_2*nu_2*(2*alpha_2**2 * K_2 + 1j*sigma)*np.exp(m_2*z)*A2
    + (1 - sigma + 2*1j*alpha_2**2)*np.exp(k*z)*C2
        
    return sigma_zz

def compute_sigma_xz_1(z,params,k,coeffs):
    
    f_ex = params[0]
    rho_1 = params[2]
    nu_1 = params[4]            
            
    omega = 2*np.pi*f_ex # in rad/s
    m_1 = np.sqrt(k**2 - 1j*omega/nu_1)
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    
    delta_1 = np.sqrt(nu_1/omega)   
    K_1 = m_1/k    
    alpha_1 = delta_1*k
    
    sigma_zz = 1j*rho_1*nu_1*alpha_1**2 * (1 + K_1**2)*np.exp(m_1*z)*A1
    + 1j*rho_1*nu_1*alpha_1**2 * (1 + K_1**2)*np.exp(-m_1*z)*B1
    + 2*alpha_1**2*np.exp(k*z)*C1
    - 2*alpha_1**2*np.exp(-k*z)*D1
        
    return sigma_zz

def compute_sigma_xz_2(z,params,k,coeffs):
    
    f_ex = params[0]
    rho_2 = params[3]
    nu_2 = params[5]            
            
    omega = 2*np.pi*f_ex # in rad/s
    m_2 = np.sqrt(k**2 - 1j*omega/nu_2)
    
    [A1,B1,C1,D1,A2,C2] = coeffs
    
    delta_2 = np.sqrt(nu_2/omega)   
    K_2 = m_2/k    
    alpha_2 = delta_2*k
    
    sigma_xz = 1j*rho_2*nu_2*alpha_2**2*(1 + K_2**2)*np.exp(m_2*z)*A2
    + 2*alpha_2**2*np.exp(k*z)*C2
        
    return sigma_xz


#%% Set parameters

h = 1e-2 # thickness of frasil in meters
rho_1 = 1000 # in kg/m3
rho_2 = 1000 # in kg/m3
nu_1 = 1e-6 # in m2/s
nu_2 = 1e-6 # in m2/s
f_ex = 3.0 # in Hz

r = rho_1/rho_2

omega = 2*np.pi*f_ex
k_linear = omega**2/g # wave vector for non dissipative linear waves 

xi_0 = 1.5e-3 # wave amplitude in meter

params = np.array([f_ex,h,rho_1,rho_2,nu_1,nu_2])

# Order of magnitudes
U_0 = xi_0*omega
p_0 = rho_1*g*xi_0
vort_0 = U_0/h

suffixe = f'h_{h:.3f}_r_{r:.1f}_nu1_{nu_1:.2e}_nu2_{nu_2:.2e}_fex_{f_ex:.2f}Hz'
sub_fig_folder = f'{fig_folder}/{suffixe}/'
if not os.path.isdir(sub_fig_folder):
    os.mkdir(sub_fig_folder)


#%% Compute determinant for different values of k 
x = np.linspace(1e-1,200,100)
y = np.linspace(1e-1,200,100)
X, Y = np.meshgrid(x, y)
Z = X + 1j * Y

# Evaluate |det(M(k))| over the grid

det_values = np.zeros_like(Z, dtype=complex)

for i in range(Z.shape[0]):
    for j in range(Z.shape[1]):
        k = Z[i, j]

        try:
            det_values[i, j] = det_M(k,params)
        except np.linalg.LinAlgError:
            det_values[i, j] = np.nan  # handle singularities or bad evals

#%% Plots determinant of matrix 

# plot absolute value

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
extents = np.array((x.min(),x.max(),y.min(),y.max()))/k_linear
c = ax.imshow(np.log10(abs(det_values)),cmap = parula_map,origin = 'lower',aspect = 'equal',
              interpolation = 'gaussian',extent = extents)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$\mathrm{log10(|det(M(k))|)}$',labelpad = 5)
ax.set_xlabel(r'$\mathrm{Re}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$\mathrm{Im}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)

ax.set_title('abs(det(M(k))) in complex plane')
fig_label = f'det_M_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Perform Newton algorithm 

###################### Parameters for Newton method #####################################
k0 = k_linear*(1 + 1j*0.1)  # initial wavevector to start Newton algorithm 
print(f'k = omega2/g leads to : {k_linear}')
# compute derivative 
f = det_M(k0,params)
df = d_det_M(k0,params,h = 1e0)
print(f/df)

tol = 1e-2
max_iter = 100
step_derivative = 1e0
threshold_derivative = 1e-80
new_k = newton_raphson(params, k0,tol,max_iter,step_derivative,threshold_derivative)
print(f'Wavevector from newton method : {new_k}')
#########################################################################################


# Plot results from Newton method
set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
extents = np.array((x.min(),x.max(),y.min(),y.max()))/k_linear
c = ax.imshow(np.log10(abs(det_values)),cmap = parula_map,origin = 'lower',aspect = 'equal',
              interpolation = 'gaussian',
              extent = extents)
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)

ax.plot(np.real(k0)/k_linear,np.imag(k0)/k_linear,'ro')
ax.plot(np.real(new_k)/k_linear,np.imag(new_k)/k_linear,'d',color = 'k')

cbar.set_label(r'$\mathrm{log10(|det(M(k))|)}$',labelpad = 5)
ax.set_xlabel(r'$\mathrm{Re}(k)/k_0 \; \mathrm{(rad.m^{-1})}$', labelpad = 5)
ax.set_ylabel(r'$\mathrm{Im}(k)/k_0\; \mathrm{(rad.m^{-1})}$', labelpad = 5)

ax.set_title('$\omega^2/g = ' + f'{k_linear:.2f}' + '\; \mathrm{rad.m^{-1}}$')

fig_label = f'Newton_method_results_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Solve system using wave amplitude as an input  

coeffs = get_coeffs(xi_0,new_k,params)

#%% Set coordinates 
############## Coordinates ####################
H = 5*h
dx = -1e-2*h
z1 = np.arange(0,-h,step = dx)
z2 = np.arange(-h,-H,step = dx)
z = np.concatenate((z1,z2))
# z1 = np.linspace(0,-h*0.99,100) # z-coordinates in meters 
# z2 = np.linspace(-h*1.01,-H,100)

# z = np.linspace(0,-H,1000)

Lx = 1.0
x = np.linspace(0,Lx,100) # x-coordinate in meter
###############################################

#%% Compute pressure 
p_1 = compute_p1(z1, new_k, coeffs)
p_2 = compute_p2(z2, new_k, coeffs)
p = np.concatenate((p_1,p_2))

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(z/h,np.real(p)/p_0,label = r'$Re(\hat{p})$')
ax.plot(z/h,np.imag(p)/p_0,label = r'$Im(\hat{p})$')

ax.set_xlabel(r'$z/h$')
ax.set_ylabel(r' $p/p_0$')
ax.legend()
fig_label = f'p_VS_z_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute vorticity fields

vort1 = compute_vort1(z1, params, new_k, coeffs)
vort2 = compute_vort2(z2, params, new_k, coeffs)
vort = np.concatenate((vort1,vort2))

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(z/h,np.real(vort)/vort_0,label = r'$Re(\hat{\Omega})$')
ax.plot(z/h,np.imag(vort)/vort_0,label = r'$Im(\hat{\Omega})$')

ax.set_xlabel(r'$z/h$')
ax.set_ylabel(r' $\Omega/\Omega_0$')
ax.legend()
fig_label = f'vort_VS_z_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute vertical profile of velocity field 

# coeffs = [A1,B1,C1,D1,A2,C2]

u1 = compute_u1(z1, params, new_k, coeffs)
u2 = compute_u2(z2, params, new_k, coeffs)

w1 = compute_w1(z1, params, new_k, coeffs)
w2 = compute_w2(z2, params, new_k, coeffs)

u = np.concatenate((u1,u2))
w = np.concatenate((w1,w2))

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(z/h,np.real(u)/U_0,label = r'$Re(\hat{u})(z)$')
ax.plot(z/h,np.imag(u)/U_0,label = r'$Im(\hat{u})(z)$')
ax.plot(z/h,np.real(w)/U_0,label = r'$Re(\hat{w})(z)$')
ax.plot(z/h,np.imag(w)/U_0,label = r'$Im(\hat{w})(z)$')

ax.set_xlabel(r'$z/h$')
ax.set_ylabel(r'$U/U_0$')
ax.legend()
fig_label = f'velocities_hat_VS_z_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

# fig_label = f'vertical_profiles_results_h_{h:.3f}_nu1_{nu_1:.2e}_nu2_{nu_2:.2e}_fex_{f_ex:.2f}_xi_{xi_0:.2e}'
# figname = f'{fig_folder}{fig_label}'
# plt.savefig(figname + '.pdf', bbox_inches='tight')
# plt.savefig(figname + '.png', bbox_inches='tight')


#%% Compute stresses

sigma_zz_1 = compute_sigma_zz_1(z1, params, new_k, coeffs)
sigma_zz_2 = compute_sigma_zz_2(z2, params, new_k, coeffs)
sigma_zz = np.concatenate((sigma_zz_1,sigma_zz_2))

sigma_xz_1 = compute_sigma_xz_1(z1, params, new_k, coeffs)
sigma_xz_2 = compute_sigma_xz_2(z2, params, new_k, coeffs)
sigma_xz = np.concatenate((sigma_xz_1,sigma_xz_2))

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(z/h,np.real(sigma_zz)/p_0,label = r'$Re(\hat{\sigma}_{zz})(z)$')
ax.plot(z/h,np.imag(sigma_zz)/p_0,label = r'$Im(\hat{\sigma}_{zz})(z)$')
ax.plot(z/h,np.real(sigma_xz)/p_0,label = r'$Re(\hat{\sigma}_{xz})(z)$')
ax.plot(z/h,np.imag(sigma_xz)/p_0,label = r'$Im(\hat{\sigma}_{xz})(z)$')

ax.set_xlabel(r'$z/h$')
ax.set_ylabel(r'$\sigma/p_0$')
ax.legend()
fig_label = f'stresses_hat_VS_z_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

#%% Compute velocity 2D velocity field

U = np.outer(u,np.exp(1j*new_k*x))
W = np.outer(w,np.exp(1j*new_k*x))

extents = np.array([x.min()/Lx,x.max()/Lx,z.min()/h,z.max()/h])

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.imshow(np.real(U)/U_0,cmap = parula_map, interpolation = 'gaussian',
          aspect = 'auto',extent = extents)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$U(x,z)/U_0$',labelpad = 5)
ax.set_xlabel(r'$x/L$')
ax.set_ylabel(r'$z/h$')
fig_label = f'U_VS_xz_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.imshow(np.real(W)/U_0,cmap = parula_map, interpolation = 'gaussian',
          aspect = 'auto',extent = extents)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$W(x,z)/U_0$',labelpad = 5)
ax.set_xlabel(r'$x/L$')
ax.set_ylabel(r'$z/h$')
fig_label = f'W_VS_xz_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


#%% Show horizontal velocity profile 

# get index in the middle of fluid 1
idx_1 = np.argmin(abs(z + h))

idx_2 = np.argmin(abs(z + 2*h))

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(x/Lx,np.real(U[idx_1,:])/U_0,label = r'$Re(u)(x,z = -h/2)$')
ax.plot(x/Lx,np.real(U[idx_2,:])/U_0,label = r'$Re(u)(x,z = -2h)$')

ax.set_xlabel(r'$x/L$')
ax.set_ylabel(r'$U/U_0$')
# ax.set_ylim([-0.8,0.8])
ax.legend()
fig_label = f'U_VS_x_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
ax.plot(x/Lx,np.real(W[idx_1,:])/U_0,label = r'$Re(w)(x,z = -h/2)$')
ax.plot(x/Lx,np.real(W[idx_2,:])/U_0,label = r'$Re(w)(x,z = -2h)$')

ax.set_xlabel(r'$x/L$')
ax.set_ylabel(r'$W/U_0$')
# ax.set_ylim([-0.8,0.8])
ax.legend()
fig_label = f'W_VS_x_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


#%% Compute stresses in 2D 

Sigma_zz = np.outer(sigma_zz,np.exp(1j*new_k*x))
Sigma_xz = np.outer(sigma_xz,np.exp(1j*new_k*x))

extents = np.array([x.min()/Lx,x.max()/Lx,z.min()/h,z.max()/h])

set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.imshow(np.real(Sigma_zz)/p_0,cmap = parula_map, interpolation = 'gaussian',
          aspect = 'auto',extent = extents)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$\sigma_{zz}(x,z)/p_0$',labelpad = 5)
ax.set_xlabel(r'$x/L$')
ax.set_ylabel(r'$z/h$')
fig_label = f'Sigma_zz_VS_xz_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')


set_graphs.set_matplotlib_param('single')
fig, ax = plt.subplots()
c = ax.imshow(np.real(Sigma_xz)/p_0,cmap = parula_map, interpolation = 'gaussian',
          aspect = 'auto',extent = extents)

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="2%", pad=0.1)

cbar = plt.colorbar(c,cax = cax)
cbar.set_label(r'$\sigma_{xz}(x,z)/p_0$',labelpad = 5)
ax.set_xlabel(r'$x/L$')
ax.set_ylabel(r'$z/h$')
fig_label = f'Sigma_xz_VS_xz_{suffixe}'
figname = f'{sub_fig_folder}{fig_label}'
plt.savefig(figname + '.pdf', bbox_inches='tight')
plt.savefig(figname + '.png', bbox_inches='tight')













