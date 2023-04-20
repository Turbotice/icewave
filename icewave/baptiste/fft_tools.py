# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 15:29:58 2023

@author: Banquise
"""
import numpy as np


def add_padding (data, padding) :
    # padding = [2**9,2**9,2**11]
    datasize = data.shape
    data_pad = np.zeros(np.power(2,padding))
    if len(padding) == 1 :
        data_pad[:datasize] = data
    if len(padding)  == 2 :
        data_pad[:datasize[0], :datasize[1]] = data
    if len(padding)  == 3 :
        data_pad[:datasize[0], :datasize[1], :datasize[2]] = data
    return data_pad