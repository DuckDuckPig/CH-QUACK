#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 14:14:42 2024

@author: jgra
"""

# In[1]:
# Import Libraries and Tools
import os
import sys
import pandas as pd
import numpy as np
import skimage.morphology
import skimage.measure
from astropy.io import fits
import sunpy.map
from reproject import reproject_interp

import warnings
warnings.filterwarnings("ignore")

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweSaveSeg_v5, acweRestoreScale

# In[2]:
# Key Variables

# Dataset
# Dataset Folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/'
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2133'

# SaveFolder
GeneralFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/EUV_SingleChannel/'
segFolder     = GeneralFolder + 'Standard/' # 'Scaled/'
uniDataFolder = GeneralFolder + 'StandardStats/' # 'ScaledStats/'

# ACWE Choice
acwePrefix  = 'ACWE.' # 'ACWEresize_param1.' # Prefix for ACWE
acweChoice  = '193'
magnetogram = 'Magnetogram'

# Save Result to
uniPrefix = 'Unipolarity.'
overwrite = False    # Run from begining 

# Inform User
verbose = True

# In[3]:
# Open File and Get List of Images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image Group we Want 
acweChoice  = np.where(keys == acweChoice)[0][0]
magnetogram = np.where(keys == magnetogram)[0][0]

# In[4]:
# Prepare to Calculate Skewness

# Find Folders
crSaveFolder   = segFolder     + CR + '/'
crSourceFolder = dataFolder    + CR + '/'
crUniFolder    = uniDataFolder + CR + '/'

# Create any Missing Folders
if not os.path.exists(crUniFolder):
    os.makedirs(crUniFolder)

# In[5]:
# Calculate Skewness

# Cycle Through Dataset
for i in range(len(data)):#-1,-1,-1):
    
    # Inform User
    if verbose:
        print(data[keys[0]][i])
    
    # Get AIA File Name
    file = data[keys[acweChoice]][i]
    
    # Placement Folders
    acweFolder = file.split('/')[0] + '/'
    if not os.path.exists(crUniFolder + acweFolder):
        os.mkdir(crUniFolder + acweFolder)
    
    # ACWE Name
    acweFile = acwePrefix + os.path.basename(file)
    acweFile = acweFolder + acweFile + '.npz'
    
    # Result Name
    uniFile = uniPrefix  + acwePrefix + os.path.basename(file)
    uniFile = acweFolder + uniFile + '.npz'
    
    # Check for File
    if overwrite or not os.path.exists(crUniFolder + uniFile):
        
        # Inform User
        if verbose:
            print('    Opening and resizing',os.path.basename(acweFile))
        
        # Open Confidence Map
        H,AH,SEG = acweSaveSeg_v5.openSeg(crSaveFolder+acweFile)
        
        # Upscale
        ResizeParam = AH['RESIZE_PARAM']
        if ResizeParam != 1:
            SEG = acweRestoreScale.upscale(SEG, AH)
        
        # # Combine
        # SEG[np.where(np.isnan(SEG))]=0
        # SEG = np.sum(SEG,axis=0)
        
        # Inform User
        if verbose:
            print('    Identifying Clusters')
        
        # Find and Label CH Regions
        ACWE_clusters = SEG * 1
        ACWE_clusters[np.where(ACWE_clusters>0)] = 1
        ACWE_clusters = skimage.morphology.dilation(ACWE_clusters.astype(bool),
                                                    np.ones([40,40]))
        ACWE_clusters = skimage.measure.label(ACWE_clusters,connectivity=2)
        
        # Split CH Regions
        SEG2 = np.zeros([np.max(ACWE_clusters),len(SEG),len(SEG[0])])
        for j in range(1,np.max(ACWE_clusters)+1):
            cluster = ACWE_clusters * 1
            cluster[np.where(cluster!=j)] = 0
            cluster = cluster/np.max(cluster)
            SEG2[j-1] = SEG*cluster
        
        # Inform User
        if verbose:
            print('    Preparing Magnetogram')
        
        # Open Magnetogram
        hdulist = fits.open(crSourceFolder + data[keys[magnetogram]][i])
        hdulist.verify('silentfix') #necessary for successful data read
        h_mag = hdulist[1].header
        J_mag = hdulist[1].data
        hdulist.close()
        
        # Create Maps
        mapHMI = sunpy.map.Map((J_mag,h_mag))
        mapHMI.plot_settings['cmap']='hmimag'
        mapACWE = sunpy.map.Map(SEG,H)
        
        # Reproject Magnetogram
        hmiReproject,footprint = reproject_interp(mapHMI,mapACWE.wcs,
                                                  mapACWE.data.shape)
        
        # Generate Weights to Address Projection Effects
        # Based on: https://docs.sunpy.org/en/stable/generated/gallery/map_transformations/reprojection_aia_euvi_mosaic.html#improving-the-output
        Weights = sunpy.map.all_coordinates_from_map(mapACWE)
        Weights = Weights.transform_to('heliocentric').z.value
        Weights = (Weights/np.nanmax(Weights)) ** 3
    
        # Make Weighted Magentogram Map
        hmiReprojectWeighted = hmiReproject * Weights
        
        # Inform User
        if verbose:
            print('    Calculating Unipolarity')
        
        # Prepare to Calculate Unipolarity
        uni = np.empty([len(SEG2),2])
        uni[:] = np.nan
        
        # Cycle Through Regions
        for j in range(len(SEG2)):
                
            # Threshold by Confidence
            CH_kept = SEG2[j]*1
            CH_kept = CH_kept.astype(bool)
            
            # Flatten Selected Region
            flattenNormal = hmiReproject[CH_kept].flatten()
            flattenWeight = hmiReprojectWeighted[CH_kept].flatten()
            
            # Remove NAN Values, if any Persist
            flattenNormal = flattenNormal[np.logical_not(np.isnan(flattenNormal))]
            flattenWeight = flattenWeight[np.logical_not(np.isnan(flattenWeight))]
            
            # Calculate Unipolarity
            numM = np.mean(np.abs(flattenNormal)) - np.abs(np.mean(flattenNormal))
            denM = np.mean(np.abs(flattenNormal))
            uniM = numM/denM
            numW = np.mean(np.abs(flattenWeight)) - np.abs(np.mean(flattenWeight))
            denW = np.mean(np.abs(flattenWeight))
            uniW = numW/denW
            
            # Store Results
            uni[j,0] = uniM # Unweighted Skewness
            uni[j,1] = uniW # Weighted Skewness
        
        # Inform User
        if verbose:
            print('    Saving Results')
        
        # Save Results
        np.savez_compressed(crUniFolder+uniFile,SEG2,uni)
    
# In[6]:
# End Process

print('**Process Complete**')
