#%%
import os
import pickle
import sys
import numpy as np
import matplotlib.pyplot as plt

#%%
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '../../../../..')))

#from icewave.field.drone import *
from icewave.drone.drone_projection import *
from icewave.drone.get_frame_from_video import *
#from icewave.drone.image_superpositions import *

from drone_modifvasco import *

# %% example of projection with a simple photo
disk = 'H:'

dict_alldata = {}
img_path = f'{disk}/data/0213/Drones/Bernache/All/DJI_20260213160139_0065_D.JPG'

# first step : find metadata to have angle, azimut, height etc.
get_records('0213')
# %%
