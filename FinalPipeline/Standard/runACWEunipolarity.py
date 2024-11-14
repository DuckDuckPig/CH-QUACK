#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 14:55:33 2024

@author: jgra
"""
# In[1]
# Import Libraries and Tools
import os
import sys
import numpy as np
import pandas as pd

# Root directory of the project
ROOT_DIR = os.path.abspath("../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweUniFunctions_v1 as au1, acweSaveSeg_v5 as as5
from SDO_tools import read_fits

# In[2]
# Key Varibles

# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/' # Update to reflect Dataset Location
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2133' # Update to reflect chosen Carrington Rotation

# SaveFolder - Update to refect location where data will be saved 
saveFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations'
saveFolder = saveFolder + '/FinalPipeline/Standard/Uni/'

# ACWE Prefix
acwePrefix = 'ACWE.' # Prefix for ACWE
overwrite = False    # If True run from begining 

# ACWE Parameters 
acweChoices = ['193','Magnetogram']
resize_param = 8
foreground_weight = np.array([1,0])
background_weight = np.array([1/50,0])
unipolarity_foreground_weight = np.array([0,1])
unipolarity_background_weight = np.array([0,1/50])
prefilter_foreground_weight = foreground_weight * 1
prefilter_background_weight = background_weight * 1
prefilter_unipolarity_foreground_weight = np.zeros(len(unipolarity_foreground_weight))
prefilter_unipolarity_background_weight = np.zeros(len(unipolarity_background_weight))
alpha = [0.3,-1]
narrowband = 2
N = 10
switch = 50
strell = 5
uniCap = 0.8

acweVerbose = False            # Show process mask generation 
correctLimbBrightening = True  # Correct for Limb Brightening
rollingAlpha = 0.01            # Incrementally Increase Alpha when Failed Threshold
prefilter_fillInitHoles = True # Fill holes in mask before running ACWE
fillInitHoles = False          # but not in pre-filterred mask
prefilterSeed = True           # Prefilter mask to remove non-unipolar regions

# Inform user
verbose = True  # Inform user about which image is being processed

# # Time ACWE
# timeFile = os.path.join(ROOT_DIR,'HMI_Experiments/TestSeedingMethods/Standard/')
# timeFile = timeFile + CR + '_timeMix_Uni.csv'

# In[4]
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image group we want
acweChoices = sorted(acweChoices)
for i in range(len(acweChoices)):
    acweChoices[i] = np.where(keys == acweChoices[i])[0][0]

# In[5]
# Prepare for ACWE

# Create or find save location
crSaveFolder = saveFolder + CR + '/'
if not os.path.exists(crSaveFolder):
    os.makedirs(crSaveFolder)

# # Prepare time file
# if not os.path.exists(timeFile):
#     with open(timeFile,'w+') as f:
#         f.write('file,time\n')

# In[6]
# Perform ACWE

# Cycle Through Dataset
for i in range(len(data[keys[acweChoices[0]]])): #range(len(data[keys[acweChoices[0]]])-1, 0,-1):
    
    files = []
    
    for j in range(len(acweChoices)):
        files.append(data[keys[acweChoices[j]]][i])
    
    # Placement Folder
    acweFolder = files[0].split('/')[0] + '/'
    if not os.path.exists(crSaveFolder + acweFolder):
        os.mkdir(crSaveFolder + acweFolder)
    
    # ACWE Name
    acweFile = acwePrefix + os.path.basename(files[0])
    for j in range(1,len(acweChoices)):
        acweFile = acweFile + '.' + keys[acweChoices[j]]
    acweFile = acweFolder + acweFile + '.npz'
    
    # Check for File
    if overwrite or not os.path.exists(crSaveFolder + acweFile):
        
        # Inform User
        if verbose:
            print('Generating', os.path.basename(acweFile))
            print('    Opening Images')
        
        # Open First file
        Itmp,_,Htmp = read_fits.openAIA(dataFolder+str(CR) +'/'+files[0])
        Ishape   = np.hstack([Itmp.shape,len(files)])
        I = np.zeros(Ishape); I[:,:,0] = Itmp
        H = [];               H.append(Htmp)
        
        # Open Rest
        for j in range(1,len(files)):
            # Magnetogram
            if 'hmi.' in files[j]:
                Itmp,Htmp = read_fits.openHMI(dataFolder+str(CR) +'/'+files[j],
                                              I[:,:,0], H[0])
                Itmp[np.isnan(Itmp)] = np.nanmin(Itmp)
            # Remaining EUVs
            else:
                Itmp,_,H = read_fits.openAIA(dataFolder+str(CR) +'/'+files[j])
            # Save
            I[:,:,j] = Itmp
            H.append(Htmp)
        
        # Inform user
        if verbose:
            print('    Running ACWE')
            
        # # Time
        # start = time.time()
        
        # Run ACWE
        seg,alphar,m = au1.run_acwe(I,H,files,resize_param,foreground_weight,
                                    background_weight,
                                    unipolarity_foreground_weight,
                                    unipolarity_background_weight,
                                    alpha,narrowband,N,
                                    prefilter_foreground_weight,
                                    prefilter_background_weight,
                                    acweVerbose,
                                    prefilter_unipolarity_foreground_weight,
                                    prefilter_unipolarity_background_weight,
                                    switch, strell,uniCap, 
                                    prefilter_fillInitHoles, 
                                    correctLimbBrightening,rollingAlpha,
                                    fillInitHoles,prefilterSeed)
        
        
        # # Time
        # end = time.time()
        # timeTotal = end-start
        # timeTotal = str(timeTotal)
        
        # Inform User
        if verbose:
            print('    Saving Results\n\n')
        
        # Save Result
        init_mask_method = 'alpha*mean(qs)'
        as5.saveSeg(crSaveFolder + acweFile,seg,H,correctLimbBrightening,resize_param,
                    [[foreground_weight,unipolarity_foreground_weight],
                     [prefilter_foreground_weight,prefilter_unipolarity_foreground_weight]],
                    [[background_weight,unipolarity_background_weight],
                     [prefilter_background_weight,prefilter_unipolarity_background_weight]],
                    m,init_mask_method,fillInitHoles,alpha,alphar,narrowband,N)
        
        # # Time
        # row = acweFile + ',' + timeTotal + '\n'
        # with open(timeFile,'a+') as f:
        #     f.write(row)
        
# In[7]
# End Process

print('**Process Complete**')