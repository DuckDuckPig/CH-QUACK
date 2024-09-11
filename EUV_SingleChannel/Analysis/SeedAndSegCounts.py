#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 12:13:34 2024

@author: jgra
"""

# In[1]
# Import Libraries and Tools
import os
import sys
import glob
import numpy as np
import skimage.morphology

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweSaveSeg_v5

# In[2]
# Key Varibles

# Dataset
segFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/'
segFolder = segFolder + 'EUV_SingleChannel/Scaled/'
CR = 'CR2133'

# ACWE Choice
acwePrefix  = 'ACWEresize_param1.' # Prefix for ACWE
acweChoice  = '193'

# Program Output
seedPrefix = 'Seed'
segPrefix  = 'Seg'
sufix      = 'Counts'

# Evaluation
strell = 40 # Size of Structring element, must be uniform across all programs for valid results

# Inform User
verbose = True

# In[3]
# Get Dataset

# File List
files = sorted(glob.glob(segFolder + CR + '/*/' + acwePrefix + '*.' + acweChoice + '.*.npz'))

# In[4]
# Prepare for Loop

# Placeholder
seedCounts = np.empty(len(files)); seedCounts[:] = np.nan
segCounts  = np.empty(len(files)); segCounts[:]  = np.nan

# In[5]
# Cycle through Dataset

for i in range(len(files)):
    
    # Get File name
    file = files[i]
    
    # Inform User
    if verbose:
        print(os.path.basename(file))
    
    # Open File
    H,AH,SEG = acweSaveSeg_v5.openSeg(file)
    initMask = AH['INIT_MASK']
    ResizeParam = AH['RESIZE_PARAM']
    
    # Inform User
    if verbose:
        print('    Identifying Clusters In Seed')
    
    # Find and Locate Seed Clusters
    SeedClusters = initMask * 1
    SeedClusters = skimage.morphology.dilation(SeedClusters.astype(bool),
                                               np.ones([strell//ResizeParam,
                                                        strell//ResizeParam]))
    SeedClusters = skimage.measure.label(SeedClusters,connectivity=2)
    seedCounts[i] = np.max(SeedClusters)
    
    # Inform User
    if verbose:
        print('    Identifying Clusters In Segmentations')
    
    # Find and Locate Set Clusters
    SegClusters = SEG * 1
    SegClusters = skimage.morphology.dilation(SegClusters.astype(bool),
                                              np.ones([strell//ResizeParam,
                                                       strell//ResizeParam]))
    SegClusters = skimage.measure.label(SegClusters,connectivity=2)
    segCounts[i] = np.max(SegClusters)
    
    if verbose:
        print()
    
# In[6]
# Save Results

# Inform Users
if verbose:
    print('Saving Results')

# Save
if not os.path.exists('Results/'):
    os.mkdir('Results')
np.savez_compressed('Results/' + seedPrefix + str(ResizeParam) + sufix + '.' + CR + '.npz', seedCounts)
np.savez_compressed('Results/' + segPrefix  + str(ResizeParam) + sufix + '.' + CR + '.npz', segCounts)

if verbose:
    print()

# In[7]
# End Process

print('**Process Complete**')
