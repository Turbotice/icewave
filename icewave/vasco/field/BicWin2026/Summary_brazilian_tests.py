#%%
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

#%%
disk = "H"
date = '0223'
# Load the data
data_path = f'{disk}:/data/{date}/Tests_Bresiliens/Summary_samples_tests.csv'

df = pd.read_csv(data_path, sep=';', encoding='latin-1')
#%%
# Display the first few rows of the DataFrame
print(df.head())
# %%
