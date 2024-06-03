#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 14:43:39 2024

@author: jgra
"""

# In[1]:
# Import Libraries and Tools
import os
import sys
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import scipy as sp
from astropy.io import fits
import sunpy.map
from aiapy.calibrate import register, update_pointing
import warnings
warnings.filterwarnings("ignore")

# Root directory of the project
ROOT_DIR = os.path.abspath("../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweFunctions_v6, acweSaveSeg_v5
from ACWE_python_fall_2023.ACWE_python_v3 import correct_limb_brightening as clb, acwe

# import time

# In[2]:
# Key Variables

# Dataset
# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/' # Update to reflect Dataset Location
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2133' # Update to reflect chosen Carrington Rotation

# SaveFolder - Update to refect location where data will be saved 
saveFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/EUV_SingleChannel/Standard_history/'

# ACWE Prefix
acwePrefix = 'ACWE.' # Prefix for ACWE
overwrite = False    # If True run from begining 

# ACWE Parameters 
acweChoice = '193'
noScale = True
resize_param = 8          # These values are the default values taken from:
foreground_weight = 1     # L. E. Boucheron, M. Valluri, and R. T. J. McAteer, 
background_weight = 1/50. # "Segementation of Coronal Holes Using Active 
alpha = 0.3               # Contours Without Edges," Solar Physics, 
narrowband = 2            # vol. 291, pp. 2353-2372, 2016
N = 10

acweVerbose = False           # Show process mask generation 
correctLimbBrightening = True # Correct for Limb Brightening
rollingAlpha = 0.01           # Incrementally Increase Alpha when Failed Threshold
fillInitHoles=True            # Fill holes in mask before running ACWE

# Inform user
verbose = True  # Inform user about which image is being processed

# ACWE on 211 Data
# acweChoice = 211; alpha = 0.3; background_weight = 1/100. # Old - refine parameters

# # Time ACWE
# timeFile = os.path.join(ROOT_DIR,'Standard/') + CR + '_timeStandard.csv'

# In[3]
# Key Functions

def itterate_acwe(I,im_size,sd_mask,m,foreground_weight=1,
                  background_weight=1/50.,narrowband=2,N=10,
                  fillInitHoles=True,verbose=False):
    
    if verbose:
        plt.ion() # interactive plotting on 
    
    # Set up variables for ACWE iterations
    if fillInitHoles:
        m_seg = sp.ndimage.morphology.binary_fill_holes(m) # fill holes
    else:
        m_seg = m # image to keep track of current initialization of ACWE
        
    history = m_seg * 1
        
    # Valid Mask - Perform ACWE
    if np.sum(m.astype(int))!=0:
        
        I_seg = I # image of current segmentation
        I_seg[~sd_mask] = I[~m_seg&sd_mask].mean() # set pixels outside SD to
                                                   # mean of background to 
                                                   # force ACWE to ignore
        counter = 0 # to keep track of proxy of iterations
        seg_diff_cum = np.zeros(im_size) # to keep track of how many times
                                         # pixels change classes over
                                         # iterations
        iterate = 1 # flag to continue iterating
        
        # Continue iterating with N iterations of the ACWE evolution, 
        # followed by check for convergence. Running ACWE for default of N=10
        # iterations is a good trade-off between checking too often and not
        # often enough for convergence.
        
        if verbose:
            print('% Diff            % New Diff')
        while iterate:
            for i in range(N):
                seg = acwe.acwe(I_seg,m_seg,1,(0,0,foreground_weight,
                                background_weight),narrowband,verbose) # Evolve ACWE 
                                                                       # for N
                                                                       # Iterations
                history = np.dstack([history,seg*1])
            
            # update current segmentation
            I_seg = I 
            I_seg[~sd_mask] = I[~seg&sd_mask].mean() 
            
            # difference in seg from previous iteration to now
            seg_diff = seg.astype(int) - m_seg.astype(int) 
            
            m_seg = seg # update m_seg image
            counter = counter + 1 # iterate counter
            
            # compute percentage of pixels that changed between previous 
            # iteration, and now
            percent_diff = float(abs(seg_diff).sum())/float(seg.sum())*100.0
            # keep track of how many times pixels have changed classes
            seg_diff_cum = seg_diff_cum + abs(seg_diff)
            # percentage of currently new pixels that have never changed
            # classes before
            percent_new_diff = float((seg_diff_cum*(abs(seg_diff)==1)).sum())/\
                               float((seg_diff_cum*(abs(seg_diff))>=1).sum()+\
                               np.finfo(float).eps)*100  
            if verbose:
                print(str(percent_diff) + ' ' + str(percent_new_diff))
            if percent_new_diff==0 | ~(seg.sum()>0):
                iterate = 0
        
        # Return Segmentation
        return seg,history
    else:
        return m * 1,history

# Run ACWE
def run_acwe(J,h,resize_param=8,foreground_weight=1,background_weight=1/50.,
             alpha=0.3,narrowband=2,N=10,verbose=False,
             correctLimbBrightening=True,rollingAlpha=0,fillInitHoles=True):
    
    # Resize image
    I,im_size,sun_radius,sun_center = acweFunctions_v6.resize_EUV(J,h,resize_param)
    
    # Correct limb brightening per Verbeeck et al. 2014
    if correctLimbBrightening:
        I = clb.correct_limb_brightening(I,sun_center,sun_radius)

    #  Define solar disk mask and initial mask
    if rollingAlpha != 0:
        sd_mask,m,alphar = acweFunctions_v6.inital_masks(I,im_size,sun_radius,
                                                         sun_center,alpha,
                                                         rollingAlpha)
    else:
        sd_mask,m = acweFunctions_v6.inital_masks(I,im_size,sun_radius,
                                                  sun_center,alpha,
                                                  rollingAlpha)
    
    # Perform ACWE
    seg,history = itterate_acwe(I,im_size,sd_mask,m,foreground_weight,
                                background_weight,narrowband,N,fillInitHoles,verbose)
    
    # Return Results
    if rollingAlpha != 0:
        return seg,alphar,history
    
    else:
        return seg,history

# In[3]:
# Open file and get list of images

# Open Carrington Rotation Document
CarringtonFile = traceFolder + CR + '.csv'
data = pd.read_csv(CarringtonFile,header=0)
keys = data.keys()

# Find Image group we want 
acweChoice = np.where(keys == acweChoice)[0][0]

# In[4]:
# Prepare for ACWE

# Create or find save location
crSaveFolder = saveFolder + CR + '/'
if not os.path.exists(crSaveFolder):
    os.makedirs(crSaveFolder)
    
# # Prepare time file
# if not os.path.exists(timeFile):
#     with open(timeFile,'w+') as f:
#         f.write('file,time\n')
    
# In[5]:
# Perform ACWE

# Cycle Through Dataset
for file in data[keys[acweChoice]]:#[len(data[keys[acweChoice]])-1:0:-1]:
    
    # Placement Folder
    acweFolder = file.split('/')[0] + '/'
    if not os.path.exists(crSaveFolder + acweFolder):
        os.mkdir(crSaveFolder + acweFolder)
    
    # ACWE Name
    acweFile = acwePrefix + os.path.basename(file)
    acweFile = acweFolder + acweFile + '.npz'
    
    # Check for File
    if overwrite or not os.path.exists(crSaveFolder + acweFile):
        
        # Inform User
        if verbose:
            print('Generating', os.path.basename(acweFile))
            print('    Opening EUV Image')
        
        # Attempt to open and update file
        success = False
        while not success:
            try:
                # Extract Image and Header Data
                hdulist = fits.open(dataFolder+str(CR) +'/'+file)
                hdulist.verify('silentfix') # no clue why this is needed for successful data read
                h = hdulist[1].header
                J = hdulist[1].data
                hdulist.close()
                
                # Update to Level 1.5 Data Product
                if h['LVL_NUM'] < 1.5:
                    m = sunpy.map.Map((J,h))    # Create Sunpy Map
                    m = update_pointing(m)      # Update Header based on Latest Information
                    m_registrered = register(m) # Recenter and rotate to Solar North
                    I = m_registrered.data
                    # Undo Keword Renaming
                    H = dict()
                    for k in m_registrered.meta.keys(): 
                        H[k.upper()] = m_registrered.meta[k] 
                # Skip if already Level 1.5
                else:
                    # Convert header to dictionary
                    m = sunpy.map.Map((J,h)) # Create Map
                    H = dict()
                    for k in m.meta.keys():
                        h[k.upper()] = m.meta[k]
                    I = J*1 # Copy image
                success = True
            except:
                pass
        
        # Inform user
        if verbose:
            print('    Running ACWE')
            
        # # Time
        # start = time.time()
        
        # Run ACWE
        seg,alphar,m = run_acwe(I,H,resize_param,
                                foreground_weight,
                                background_weight,
                                alpha,narrowband,
                                N,acweVerbose,
                                correctLimbBrightening,
                                rollingAlpha,
                                fillInitHoles)
        
        # # Time
        # end = time.time()
        # timeTotal = end-start
        # timeTotal = str(timeTotal)
        
        # Inform User
        if verbose:
            print('    Saving Results\n\n')
        
        # Save Result
        init_mask_method = 'alpha*mean(qs)'
        acweSaveSeg_v5.saveSeg(crSaveFolder + acweFile,seg,H,
                               correctLimbBrightening,resize_param,
                               foreground_weight,background_weight,m,
                               init_mask_method,fillInitHoles,alpha,
                               alphar,narrowband,N)
        
        # # Time
        # row = acweFile + ',' + timeTotal + '\n'
        # with open(timeFile,'a+') as f:
        #     f.write(row)
        
# In[6]:
# End Process

print('**Process Complete**')
