#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 13:37:08 2024

@author: jgra
"""

# In[1]
# Import Libraries and Tools
import os
import sys
import pandas as pd
import numpy as np
import glob
import skimage.morphology

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweSaveSeg_v5#, acweRestoreScale

# In[2]
# Key Varibles

# Dataset
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR = 'CR2101'  # Update to reflect chosen Carrington Rotation

# ACWE Folders
ACWEmixFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/HMI_Experiments/TestSeedingMethods/Standard/Uni/'
compareFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/HMI_Experiments/TestEvolutionMethods/Standard/Uni/'

# ACWE Parameters 
acweChoice = '193'
hmiChoice  = 'Magnetogram'

# Inform user
verbose = True

# Tracibility
fileSufix = '.HMItypeDiscrepancies.csv'

# In[3]
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image groups we want 
acweChoice = np.where(keys == acweChoice)[0][0]
hmiData    = np.where(keys == hmiChoice)[0][0]

# In[4]
# Prepare for analysis

# Create or find save location
crSaveFolder1 = ACWEmixFolder  + CR + '/'
crSaveFolder2 = compareFolder  + CR + '/'

# Overwrite File
resultsFolder = 'HMI_Experiments/TestSeedingMethods/Analysis/Results/'
resultsFolder = os.path.join(ROOT_DIR,resultsFolder)
if not os.path.exists(resultsFolder):
    os.makedirs(resultsFolder)
traceFile = resultsFolder + CR + fileSufix
with open(traceFile,'w+') as f:
    f.write('Index,File Name,Old,New\n')


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
    _,AM,SegM = acweSaveSeg_v5.openSeg(mixImage)
    _,AC,SegC = acweSaveSeg_v5.openSeg(cmpImage)
    
    # # Resize Segmentations
    # SegM = acweRestoreScale.upscale(SegM,AM)
    # SegC = acweRestoreScale.upscale(SegC,AC)
    
    # Label CHs
    segm = skimage.morphology.dilation(SegM.astype(bool),
                                       np.ones([5,5]))#np.ones([40,40]))
    segm = skimage.measure.label(segm,connectivity=2)
    
    segc = skimage.morphology.dilation(SegC.astype(bool),
                                       np.ones([5,5]))#np.ones([40,40]))
    segc = skimage.measure.label(segc,connectivity=2)
    
    # Determin Number of Regions
    Mcount = np.max(segm)
    Ccount = np.max(segc)
    
    # Report Number of Regions
    if Mcount != Ccount:
        entry  = str(i)+',' + os.path.basename(file) + ','\
               + str(Ccount) + ',' + str(Mcount) + '\n'
        with open(traceFile,'a+') as f:
            f.write(entry)
            
# In[6]
# End Program
print('**Process Complete**')