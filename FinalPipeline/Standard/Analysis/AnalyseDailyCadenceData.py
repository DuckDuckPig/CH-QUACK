#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 10:45:58 2025

@author: jgra
"""

# In[1]
# Import Libraries and Tools

import os
import sys
import glob
import sunpy.map
import numpy as np

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweSaveSeg_v5 as as5, acweRestoreScale as ars
from DatasetTools import DataManagmentTools as dmt

# In[2]
# Key Varibles

# Datset
dataFolder  = '/home/jgra/Coronal Holes/DailyCadenceData/'
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DailyCadenceDownloadLists/')

# SaveFolder - Update to refect location where data will be saved 
saveFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations'
saveFolder = saveFolder + '/FinalPipeline/DailyCadence/Standard/Uni/'

# Save Result to
saveName = '.Area.npz'

# Inform User
verbose = True

# Overwrite existing data
overwrite = False

# In[3]
# Prepare to Cycle Through Dataset

# Get List of CRs/Datset Folders
CRsaveFolders = sorted(glob.glob(saveFolder + '*'))

resultsFolder = 'Results/'
if not os.path.exists(resultsFolder):
    os.makedirs(resultsFolder)

# In[4]
# Cycle Through Dataset

for i in range(len(CRsaveFolders)):
    
    CRsaveFolder = CRsaveFolders[i]
    
    # Get CR
    CR = os.path.basename(CRsaveFolder)
    
    # Determin File Name
    results = resultsFolder + CR + saveName
    
    # Cycle through entries
    if overwrite or not(os.path.exists(results)):
    
        # Inform user
        if verbose:
            print(CR)
        
        # Get List of images
        files = sorted(glob.glob(CRsaveFolder + '/*/*.npz'))
        
        # Prep
        dates   = []
        pxAreas = []
        MmAreas = []
        
        # Cycyle Through CR
        for file in files:
            
            # Inform User
            if verbose:
                print('   ', os.path.basename(file))
                print('        Opening File')
            
            # Open File
            H,AH,SEG = as5.openSeg(file)
            
            # Inform User
            if verbose:
                print('        Restoring Scale')
            
            # Resize File
            SEG = ars.upscale(SEG,AH)
            
            # Inform User
            if verbose:
                print('        Gathering Coordinate Data')
                
            # Convert to sunpy Map
            mapACWE = sunpy.map.Map((SEG,H[0]))
            
            # Get Coordinate Information
            cord = sunpy.map.all_coordinates_from_map(mapACWE)
            
            # Convert to heliocentric
            cord = cord.transform_to('heliocentric')
            
            # Corridante Transform (Spherical)
            theta = np.arctan((np.sqrt(cord.x**2+cord.y**2)/cord.z))
            
            # Inform User
            if verbose:
                print('        Calculating Areas')
            
            # Area in Mm
            area   = 0.189/np.cos(theta)
            MmAreas.append(np.nansum(area*SEG.astype(int)))
            
            # Areas in pixels
            pxAreas.append(np.nansum(SEG.astype(int)))
            
            # date
            dates.append(dmt.timeUnformat(H[0]['T_REC']))
            
            print()
    
        # Save Data
        np.savez_compressed(results, dates, pxAreas, MmAreas)

# In[5]
# End Process

print('**Process Complete**')