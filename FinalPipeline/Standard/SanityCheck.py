#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 15:21:06 2024

@author: jgra
"""

# In[1]
# Import libraries and tools
import os
import sys
import numpy as np
import glob

# Root directory of the project
ROOT_DIR = os.path.abspath("../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweSaveSeg_v5 as as5
from Metrics import JaccardIndexMetric as jim

# In[2]
# Datasets

mainFolder  = '/mnt/coronal_holes/Paper 2/Code 02 Observations/'
dataFolder1 = mainFolder + 'FinalPipeline/Standard/Uni/CR2133/'
dataFolder2 = mainFolder + 'HMI_Experiments/TestSeedingMethods/Standard/Uni/CR2133/' 

# In[3]
# Get List of images

dataset1 = sorted(glob.glob(dataFolder1 + '*/ACWE.*.npz'))
dataset2 = sorted(glob.glob(dataFolder2 + '*/ACWE.*.npz'))

# In[4]
# Check

if len(dataset1) != len(dataset2):
    
    print('Files Missing')
    
else:
    
    IOU = np.empty(len(dataset1));IOU[:] = 0
    
    for i in range(len(dataset1)):
        
        _,_,seg1 = as5.openSeg(dataset1[i])
        _,_,seg2 = as5.openSeg(dataset2[i])
        
        IOU[i] = jim.IOU(seg1,seg2)
        
print(np.nanmin(IOU))
print(np.nanmax(IOU))
print(np.count_nonzero(IOU-1))