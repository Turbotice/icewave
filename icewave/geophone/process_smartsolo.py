%matplotlib qt
import os
import shutil
from obspy import read
from obspy.core import UTCDateTime
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector
import mplcursors
from datetime import datetime
import numpy as np



path2data = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/2024_BICWIN/0210/Geophones'
acqu_numb = '0001'
geophones_table_path = '/Users/moreaul/Documents/Travail/Projets_Recherche/MSIM/data/geophones_table'

def read_data(path2data):
    miniseed_files = []
    # Iterate over files in the specified directory
    for filename in os.listdir(path2data):
        if filename.endswith(".miniseed"):
            file_path = os.path.join(path2data, filename)
            miniseed_files.append(file_path)
    # Read MiniSEED files
    streams = []
    for file_path in miniseed_files:
        stream = read(file_path)
        streams.append(stream)
    return streams, miniseed_files


seismic_data_streams, miniseed_files = read_data(path2data +'/' +acqu_numb)

