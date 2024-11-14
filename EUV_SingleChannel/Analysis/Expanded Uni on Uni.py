#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 08:38:21 2024

@author: jgra
"""
# In[1]
# Data Prep

# In[1a]
# Import Libraries and Tools
import os
import sys
from astropy.io import fits
import sunpy.map
from reproject import reproject_interp
import warnings
warnings.filterwarnings("ignore")
import glob
import numpy as np
import skimage.morphology

# ACWE utilities
# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweSaveSeg_v5, acweRestoreScale

# In[1b]
# Key Varibels 

# Carrington Rotation
CR = 'CR2133'

# Datasets
dataFolder    = '/home/jgra/Coronal Holes/newDataset/'
GeneralFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/'
evolveFolder  = GeneralFolder + 'EUV_SingleChannel/Standard_history/'
uniElvFolder  = GeneralFolder + 'HMI_Experiments/TestEvolutionMethods/Standard_history/Uni/'
stats1Folder  = GeneralFolder + 'EUV_SingleChannel/ScaledStats/'
stats8Folder  = GeneralFolder + 'EUV_SingleChannel/StandardStats/'
seeds1Folder  = GeneralFolder + 'EUV_SingleChannel/ScaledSeedStats/'
seeds8Folder  = GeneralFolder + 'EUV_SingleChannel/StandardSeedStats/'

# Seed Counts by file
seed8title = 'Results/Seed8Counts.' + CR + '.npz'

# Weight compare by Seed
seed8toseg8Utitle = 'Results/UniSeed8toseg8U.npz'
seed8toseg8Wtitle = 'Results/UniSeed8toseg8W.npz'

# Save Information
prefix = 'UniOnUniEvolve.'
folder = 'EvolveStats/'
overwrite = False

# Inform user
verbose = True

# In[1c]
# Ensure Folder Exits
if not os.path.exists(folder):
    os.makedirs(folder)
    
# In[1d]
# Get Data
hmi    = sorted(glob.glob(dataFolder   + CR + '/*/hmi.*.fits'))
aia    = sorted(glob.glob(dataFolder   + CR + '/*/*.193.*.fits'))
orig   = sorted(glob.glob(evolveFolder + CR + '/*/*.npz'))
uni    = sorted(glob.glob(uniElvFolder + CR + '/*/*.npz'))
stats8 = sorted(glob.glob(stats8Folder + CR + '/*/*.npz'))
seeds8 = sorted(glob.glob(seeds8Folder + CR + '/*/*.npz'))

# In[1e]
# Load Data

# Load Number of Seeds/Image
data  = np.load(seed8title, allow_pickle=True)
lst   = data.files
nbrSeed8 = data[lst[0]]

# Load Unweighted Compair 
data  = np.load(seed8toseg8Utitle, allow_pickle=True)
lst   = data.files
seed8toseg8U = data[lst[0]]

# Load Weighted Compair
data  = np.load(seed8toseg8Wtitle, allow_pickle=True)
lst   = data.files
seed8toseg8W = data[lst[0]]

# In[2]
# Prioritize

# In[2a]
# Calculate Difference
unweightedDifference = seed8toseg8U[:,0] - seed8toseg8U[:,1]
weightedDifference   = seed8toseg8W[:,0] - seed8toseg8W[:,1]

# Organize
sortUnweight = sorted(unweightedDifference)
sortWeight   = sorted(weightedDifference)

# In[2b]
# Choose sample sets

# Number of Samples
size = len(sortUnweight)

# Placeholders
indexUnotCH = np.empty(size); indexUnotCH[:] = np.nan
indexWnotCH = np.empty(size); indexWnotCH[:] = np.nan
indexUisCH  = np.empty(size); indexUisCH[:]  = np.nan
indexWisCH  = np.empty(size); indexWisCH[:]  = np.nan

# Get Samples
for i in range(size):
    
    # Not CH Cases, requested number of cases, organized from Worst to best
    indexUnotCH[i] = np.where(unweightedDifference==sortUnweight[i])[0][0]
    indexWnotCH[i] = np.where(weightedDifference  ==sortWeight[i]  )[0][0]
    
    # Is CH Cases, requested number of cases, organized from worst to best
    indexUisCH[i]  = np.where(unweightedDifference==sortUnweight[len(sortUnweight)-i-1])[0][0]
    indexWisCH[i]  = np.where(weightedDifference  ==sortWeight[len(sortUnweight)  -i-1])[0][0]

# In[2c]
# Placeholders
imageIndexUnotCH = np.empty(size); imageIndexUnotCH[:] = np.nan
imageIndexWnotCH = np.empty(size); imageIndexWnotCH[:] = np.nan
imageIndexUisCH  = np.empty(size); imageIndexUisCH[:]  = np.nan
imageIndexWisCH  = np.empty(size); imageIndexWisCH[:]  = np.nan

# Find Index of images
seedCumsum = np.cumsum(nbrSeed8)

for i in range(size):
    
    imageIndexUnotCH[i] = np.where(seedCumsum>indexUnotCH[i])[0][0]
    imageIndexWnotCH[i] = np.where(seedCumsum>indexWnotCH[i])[0][0]
    imageIndexUisCH[i]  = np.where(seedCumsum>indexUisCH[i])[0][0]
    imageIndexWisCH[i]  = np.where(seedCumsum>indexWisCH[i])[0][0]

# Convert to Integer Array
imageIndexUnotCH = imageIndexUnotCH.astype(int)
imageIndexWnotCH = imageIndexWnotCH.astype(int)
imageIndexUisCH  = imageIndexUisCH.astype(int)
imageIndexWisCH  = imageIndexWisCH.astype(int)

# In[2d]
# Reset seed Index
indexUnotCH = indexUnotCH - seedCumsum[imageIndexUnotCH - 1]
indexWnotCH = indexWnotCH - seedCumsum[imageIndexWnotCH - 1]
indexUisCH  = indexUisCH  - seedCumsum[imageIndexUisCH  - 1]
indexWisCH  = indexWisCH  - seedCumsum[imageIndexWisCH  - 1]

# Underflow Case
indexUnotCH[np.where(indexUnotCH < 0)] = indexUnotCH[np.where(indexUnotCH < 0)] + seedCumsum[-1]
indexWnotCH[np.where(indexWnotCH < 0)] = indexWnotCH[np.where(indexWnotCH < 0)] + seedCumsum[-1]
indexUisCH[np.where(indexUisCH < 0)]   = indexUisCH[np.where(indexUisCH < 0)]   + seedCumsum[-1]
indexWisCH[np.where(indexWisCH < 0)]   = indexWisCH[np.where(indexWisCH < 0)]   + seedCumsum[-1]

# Convert to Integer Array
indexUnotCH = indexUnotCH.astype(int)
indexWnotCH = indexWnotCH.astype(int)
indexUisCH  = indexUisCH.astype(int)
indexWisCH  = indexWisCH.astype(int)

# In[3]
# Magnetic field Data 

# In[3a]
# Key Functions

# Calculate Unipolarity
def unipolarity(CH_kept,Hu,Hw):
    
    # Flatten Selected Region
    flattenNormal = Hu[CH_kept].flatten()
    flattenWeight = Hw[CH_kept].flatten()

    # Remove NAN Values, if any Persist
    flattenNormal = flattenNormal[np.logical_not(np.isnan(flattenNormal))]
    flattenWeight = flattenWeight[np.logical_not(np.isnan(flattenWeight))]

    # Calculate Unipolarity
    numU = np.mean(np.abs(flattenNormal)) - np.abs(np.mean(flattenNormal))
    denU = np.mean(np.abs(flattenNormal))
    uniU = numU/denU
    numW = np.mean(np.abs(flattenWeight)) - np.abs(np.mean(flattenWeight))
    denW = np.mean(np.abs(flattenWeight))
    uniW = numW/denW
    
    # Return Results
    return uniU,uniW

# In[3b]
# Calculate 

for i in range(size-1,-1,-1):
    
    title = folder + prefix + str(imageIndexUnotCH[i]) + '.' + str(indexUnotCH[i]) + '.npz'
    
    if not os.path.exists(title) or overwrite:
    
        # Inform User
        if verbose:
            print(os.path.basename(uni[imageIndexUnotCH[i]]))
            print('Seed Group:', indexUnotCH[i] + 1)
            print('    Opening Segmentation')
        
        # Open Segmentation
        H,AH,SEG = acweSaveSeg_v5.openSeg(uni[imageIndexUnotCH[i]])
        History   = AH['INIT_MASK']
        
        # Upscale
        SEG = acweRestoreScale.upscale(SEG, AH)
        
        # Inform User
        if verbose:
            print('    Opening Magnetogram')
        
        # Open Magnetogram
        hdulist = fits.open(hmi[imageIndexUnotCH[i]])
        hdulist.verify('silentfix') #necessary for successful data read
        h_mag = hdulist[1].header
        J_mag = hdulist[1].data
        hdulist.close()
    
        # Create Maps
        mapHMI = sunpy.map.Map((J_mag,h_mag))
        mapHMI.plot_settings['cmap']='hmimag'
        mapACWE = sunpy.map.Map(SEG,dict(H[0]))
    
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
            print('    Evalulating Evolution')
        
        # Prepare for unipolarity calcuations
        uniU = np.empty(len(History[0,0])); uniU[:] = np.nan
        uniW = np.empty(len(History[0,0])); uniW[:] = np.nan
        area = np.empty(len(History[0,0])); area[:] = np.nan
        
        # cycle through history
        for j in range(len(History[0,0])):
            
            # Inform User
            if verbose:
                print('       ', 1+j, '/', len(History[0,0]))
            
            # Upscale System
            Step = History[:,:,j]
            Step = acweRestoreScale.upscale(Step, AH)
            
            # Find and Label CH Regions
            ACWE_cluster = Step > 0
            ACWE_cluster = skimage.morphology.dilation(ACWE_cluster.astype(bool),
                                                        np.ones([40,40]))
            
            # Save Target Seed
            if j == 0:
                ACWE_clusterLabels = skimage.measure.label(ACWE_cluster,connectivity=2)
                SEED = ACWE_clusterLabels == (indexUnotCH[i] + 1)
                SEED = SEED * Step
            
            # Find Region
            ACWE_cluster = skimage.morphology.reconstruction(SEED*ACWE_cluster,
                                                             ACWE_cluster,
                                                             method='dilation')
            ACWE_cluster = ACWE_cluster * Step
            
            # Calculate Area
            area[j] = np.sum(ACWE_cluster)
            
            # Calculate Unipolarity
            uniU[j],uniW[j] = unipolarity(ACWE_cluster.astype(bool),
                                          hmiReproject,hmiReprojectWeighted) 
        
        
        # Inform User
        if verbose:
            print('    Saving Results')
        
        # Save Results
        np.savez_compressed(title,uniU,uniW,area)
        
        # Spacer
        if verbose:
            print()