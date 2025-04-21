#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 17:26:04 2025

@author: jeremygrajeda
"""

# In[1]
# Import Libarary and Tools
import sys
import os
import glob
import numpy as np
from scipy.io import readsav

# Root directory of the project
ROOT_DIR = os.path.abspath("../../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweUniFunctions_v1 as au1, acweSaveSeg_v5 as as5

# In[2]
# Key Varibles

# Dataset
euvDataFolder = sorted(glob.glob(ROOT_DIR+'/Extensions/ISWAT/Dataset/Uncompressed/193/CHdataprocessedkm193*.save'))
hmiDataFolder = sorted(glob.glob(ROOT_DIR+'/Extensions/ISWAT/Dataset/Uncompressed/HMI/CHboundariesHMImagnetogramspart*.save'))

# SaveFolder - Update to refect location where data will be saved 
saveFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations'
saveFolder = saveFolder + '/FinalPipeline/ISWAT/Standard/'

# ACWE Prefix
acwePrefix = 'QUACK.' # Prefix for ACWE
overwrite = False    # If True run from begining 

# ACWE Parameters 
acweChoices = ['193.','hmi.']
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

if rollingAlpha == 0:
    saveFolder = saveFolder + 'alpha' + str(int(alpha[0]*10)) + '/'
else:
    saveFolder = saveFolder + 'Results' + '/'

# In[3]
# Get Data

# Placeholder
euv_data = []
hmi_data = []

# Loop through Data
for i in range(len(euvDataFolder)):
    
    # Grab Files
    euvFile = euvDataFolder[i]
    hmiFile = hmiDataFolder[i]
    
    # Open Files
    euv_sav_data = readsav(euvFile)
    hmi_sav_data = readsav(hmiFile)
    euv_lst = list(euv_sav_data.keys())
    hmi_lst = list(hmi_sav_data.keys())
    
    # Unpack and Combine Data
    euv_data = [*euv_data,*euv_sav_data(euv_lst[0])]
    hmi_data = [*hmi_data,*hmi_sav_data(hmi_lst[0])]

# In[4]
# Prepare for ACWE

# Create or find save location
if not os.path.exists(saveFolder):
    os.makedirs(saveFolder)

# In[5]
# Perform ACWE

# Cycle Through Dataset
for i in range(len(euv_data)):
    
    #### EUV ####
    # Get EUV Entry
    euv_datum = euv_data[i]
    
    # Seperate EUV Image and Data
    Jeuv = euv_datum[0]
    Heuv = []
    for j in range(1,len(euv_datum)):
        Heuv.append(euv_datum[j])
        
    # Jerry Rig for QUACK - EUV Data
    euv_im_dims = np.array(Jeuv.shape)
    c = (euv_im_dims/2) + 1
    euvSource = Heuv[5].decode('ascii')
    hEUV = {'RSUN':Heuv[14],
            'CDELT1':Heuv[2],
            'CDELT2':Heuv[3],
            'T_REC':Heuv[4].decode('ascii'),
            'CRPIX1':c[1],
            'CRPIX2':c[0]}
    
    #### HMI ####
    # Get HMI Entry
    hmi_datum = hmi_data[i]
    
    # Seperate HMI Image and Data
    Jhmi = hmi_datum[0]
    Hhmi = []
    for j in range(1,len(hmi_datum)):
        Hhmi.append(hmi_datum[j])
    
    # Jerry Rig for QUACK - HMI Data
    hmi_im_dims = np.array(Jhmi.shape)
    c = (hmi_im_dims/2) + 1
    hmiSource = Hhmi[5].decode('ascii')
    hHMI = {'RSUN':Hhmi[14],
            'CDELT1':Hhmi[2],
            'CDELT2':Hhmi[3],
            'T_REC':Hhmi[4].decode('ascii'),
            'CRPIX1':c[1],
            'CRPIX2':c[0]}
    
    #### Jerry Rig a name for the file ####
    acweFile = acwePrefix + euvSource + '.' + hEUV['T_REC'] \
               + '.' + hmiSource + '.' + hHMI['T_REC'] + '.npz'
    acweFile = acweFile.replace(' ','_').replace(':','')
    
    # Check for File
    if overwrite or not os.path.exists(saveFolder + acweFile):
        
        # Inform User
        if verbose:
            print('Generating', os.path.basename(acweFile))
            print('    Combining Images')
        
        # Combine Images
        I = np.dstack([Jeuv,Jhmi])
        H = {0:hEUV,
             1:hHMI}
        
        # Inform user
        if verbose:
            print('    Running ACWE')
            
        if rollingAlpha == 0:
            # Run ACWE
            seg,m = au1.run_acwe(I,H,acweChoices,resize_param,foreground_weight,
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
            as5.saveSeg(saveFolder + acweFile,seg,H,correctLimbBrightening,resize_param,
                        [[foreground_weight,unipolarity_foreground_weight],
                         [prefilter_foreground_weight,prefilter_unipolarity_foreground_weight]],
                        [[background_weight,unipolarity_background_weight],
                         [prefilter_background_weight,prefilter_unipolarity_background_weight]],
                        m,init_mask_method,fillInitHoles,alpha,alpha,narrowband,N)
        else:
            # Run ACWE
            seg,alphar,m = au1.run_acwe(I,H,acweChoices,resize_param,foreground_weight,
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
            
            # Inform User
            if verbose:
                print('    Saving Results\n\n')
            
            # Save Result
            init_mask_method = 'alpha*mean(qs)'
            as5.saveSeg(saveFolder + acweFile,seg,H,correctLimbBrightening,resize_param,
                        [[foreground_weight,unipolarity_foreground_weight],
                         [prefilter_foreground_weight,prefilter_unipolarity_foreground_weight]],
                        [[background_weight,unipolarity_background_weight],
                         [prefilter_background_weight,prefilter_unipolarity_background_weight]],
                        m,init_mask_method,fillInitHoles,alpha,alphar,narrowband,N)    
        
# In[6]
# End Process

print('**Process Complete**')
