#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from load_force_vs_frame import *
from load_matdata_JZmethod import load_displacement_px
from load_displacement_vs_frame import load_disp_vs_frame

from force_disp import load_force_displacement_curve

#%%
load_force_displacement_curve(idx=0)
# %%
