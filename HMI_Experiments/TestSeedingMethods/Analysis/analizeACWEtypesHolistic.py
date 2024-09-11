#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 10:54:37 2024

@author: jgra
"""

# In[1]
# Import Libraries and Tools
import os
import sys
import pandas as pd
import numpy as np
import glob
from skimage.metrics import structural_similarity as ssim

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweSaveSeg_v5
from Metrics import consistancyErrorMetricsIII as cem, JaccardIndexMetric as jim

# In[2]
# Key Varibles

# Dataset
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR = 'CR2133'  # Update to reflect chosen Carrington Rotation

# ACWE Folders
ACWEmixFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/HMI_Experiments/TestSeedingMethods/Standard/Uni/'
#compareFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/EUV_SingleChannel/Standard/'
compareFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/HMI_Experiments/TestEvolutionMethods/Standard/Uni/'

# Tracibility
#fileSufix = '.HMItoStandardStats.npz'
fileSufix = '.HMItypeStats.npz'

# ACWE Parameters 
acweChoice = '193'

# Inform user
verbose = True

# In[3]
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image group we want 
acweChoice = np.where(keys == acweChoice)[0][0]

# In[4]
# Prepare for analysis

# Create or find save location
crSaveFolder1 = ACWEmixFolder  + CR + '/'
crSaveFolder2 = compareFolder  + CR + '/'

# Create placeholder for stats
IOU  = np.empty(len(data)); IOU[:]  = np.nan
SSIM = np.empty(len(data)); SSIM[:] = np.nan
GCE  = np.empty(len(data)); GCE[:]  = np.nan
LCE  = np.empty(len(data)); LCE[:]  = np.nan

# In[5]
# Perform analysis

# Cycle Through Dataset
for i in range(len(data[keys[acweChoice]])):
    
    # Get file
    file = data[keys[acweChoice]][i]
    
    # Inform User
    if verbose:
        print()
        print(os.path.basename(file))
        
    # Find Images
    acweFolder = file.split('/')[0] + '/'
    mixImage = glob.glob(crSaveFolder1+acweFolder+'*'+os.path.basename(file)+'*')[0]
    cmpImage = glob.glob(crSaveFolder2+acweFolder+'*'+os.path.basename(file)+'*')[0]
    
    # Open Files
    _,_,SegM = acweSaveSeg_v5.openSeg(mixImage)
    _,_,SegC = acweSaveSeg_v5.openSeg(cmpImage)
    
    # Calculate IOU
    IOU[i] = jim.IOU(SegM,SegC,True)
    
    # SSIM
    SSIM[i] = ssim(SegM,SegC)
    
    # GCE & LCE
    GCE[i],LCE[i] = cem.CE(SegM,SegC)
    
# In[6]
# Save Results for Graphing
resultsFolder = 'HMI_Experiments/TestSeedingMethods/Analysis/Results/'
resultsFolder = os.path.join(ROOT_DIR,resultsFolder)
if not os.path.exists(resultsFolder):
    os.makedirs(resultsFolder)
filename = os.path.basename(CR) + fileSufix
filename = resultsFolder + filename
np.savez(filename,IOU,SSIM,GCE,LCE)

# In[7]
# End Program
print('**Process Complete**')

    