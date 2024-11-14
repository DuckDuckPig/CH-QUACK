#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 12:46:29 2024

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
ROOT_DIR = os.path.abspath("../../../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweSaveSeg_v5
from Metrics import consistancyErrorMetricsIII as cem, JaccardIndexMetric as jim

# In[2]
# Key Varibles

# Dataset
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2133'  # Update to reflect chosen Carrington Rotation


# ACWE Folders
SaveFolderMag = '/mnt/coronal_holes/Paper 2/Code 02 Observations/'
SaveFolderMag = SaveFolderMag + 'HMI_Experiments/TestEvolutionMethods/Scaled/'
SaveFolderEUV = '/mnt/coronal_holes/Paper 1/Code Paper I Observations/Scaled/'

# ACWE Data
magMethod  = 'Uni' # 'FluxImb' 'Homogentiy' 'Skew' 
ACWEprefix = 'ACWEresize_param1.'
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
crSaveFolder1 = SaveFolderMag  + magMethod + '/' + CR + '/'
crSaveFolder2 = SaveFolderEUV  + CR + '/'

# Create placeholder for stats
outputShape = len(data)
IOU  = np.empty(outputShape); IOU[:]  = np.nan
SSIM = np.empty(outputShape); SSIM[:] = np.nan
GCE  = np.empty(outputShape); GCE[:]  = np.nan
LCE  = np.empty(outputShape); LCE[:]  = np.nan

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
        
    # Find all Images
    acweFolder  = file.split('/')[0] + '/'
    acweMagFile = glob.glob(crSaveFolder1+acweFolder+ACWEprefix+os.path.basename(file)+'*')[0]
    acweEUVfile = glob.glob(crSaveFolder2+acweFolder+ACWEprefix+os.path.basename(file)+'*')[0]
    
    # Open Images
    _,_,SegM = acweSaveSeg_v5.openSeg(acweMagFile)
    _,_,SegE = acweSaveSeg_v5.openSeg(acweEUVfile)
    
    # Calculate IOU
    IOU[i] = jim.IOU(SegM,SegE,True)
    
    # SSIM
    SSIM[i] = ssim(SegM,SegE)
    
    # GCE & LCE
    GCE[i],LCE[i] = cem.CE(SegM,SegE)
    
# In[6]
# Save Results for Graphing
resultsFolder = 'HMI_Experiments/TestEvolutionMethods/Scaled/Analysis_EUV/Results/'
resultsFolder = os.path.join(ROOT_DIR,resultsFolder)
if not os.path.exists(resultsFolder):
    os.mkdir(resultsFolder)
filename = os.path.basename(CarringtonFile) + '.' + magMethod + '.MagToEUV.npz'
filename = resultsFolder + filename
np.savez(filename,IOU,SSIM,GCE,LCE)

# In[7]:
# End Program
print('**Process Complete**')