#%%
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.optimize import curve_fit

from load_press_data import load_single_test_curve
sys.path.append('../../../tools')
from clickonfigures import get_n_points



#%%
acqnum = '2'
%matplotlib inline
dict_single_test = load_single_test_curve(acqnum=acqnum,plot=True)

Actuator_mm = dict_single_test['Actuator'].values
Load_kN = dict_single_test['Load'].values
Time_sec = dict_single_test['Time'].values

Actuator_m = 1e-3 * Actuator_mm
Load_N = 1e3 * Load_kN


n_points = int(input('numbers of clicks ?'))
%matplotlib qt
coords = get_n_points(Actuator_mm, Load_kN, n_points=n_points)

# récuperer deulement la coord x et en déduire y par interpolation 
# de la courbe de Force sv displacement
# ensuite fitter entre les deux points avec une fonciton affine

# %%
