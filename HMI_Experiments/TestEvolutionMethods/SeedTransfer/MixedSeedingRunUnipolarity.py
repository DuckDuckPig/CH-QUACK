#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 15:46:27 2024

@author: jgra
"""

# In[1]
# Import Libraries and Tools
import os
import sys
import scipy as sp
import numpy as np
from astropy.io import fits
import sunpy.map
from aiapy.calibrate import update_pointing, register
import pandas as pd
import skimage

# import warnings
# warnings.filterwarnings("ignore")

# Root directory of the project
ROOT_DIR = os.path.abspath("../../../")

# Import ACWE Tools
sys.path.append(ROOT_DIR)
from ACWE_python_fall_2023 import acweFunctions_v6 as af6, acweSaveSeg_v5 as as5
from ACWE_python_fall_2023.ACWE_python_v3 import correct_limb_brightening as clb, acwe_exp2

# import time

# In[2]
# Key Variables

# Dataset
# Dataset folders
dataFolder  = '/home/jgra/Coronal Holes/newDataset/' # Update to reflect Dataset Location
traceFolder = os.path.join(ROOT_DIR, 'DatasetTools/DownloadLists/')
CR          = 'CR2133' # Update to reflect chosen Carrington Rotation

# SaveFolder - Update to refect location where data will be saved 
saveFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/HMI_Experiments/'
saveFolder = saveFolder + 'TestEvolutionMethods/SeedTransfer/Uni/'

# Seed Sorce
seedFolder = '/mnt/coronal_holes/Paper 2/Code 02 Observations/HMI_Experiments/'
seedFolder = seedFolder + 'TestEvolutionMethods/Scaled/Uni/'
seedPrefix = 'ACWEresize_param1.'

# ACWE Prefix
acwePrefix = 'ACWE.' # Prefix for ACWE
overwrite = False    # If True run from begining 

# ACWE Parameters 
acweChoices = ['193','Magnetogram']
noScale = True
resize_param = 8
foreground_weight = np.array([1,0])
background_weight = np.array([1/50,0])
inunipolarityweight  = np.array([0,1])
outunipolarityweight = np.array([0,1/50])
#alpha = [0.3,-1]
narrowband = 2
N = 10

acweVerbose = False           # Show process mask generation 
correctLimbBrightening = True # Correct for Limb Brightening
rollingAlpha = 0.01           # Incrementally Increase Alpha when Failed Threshold
fillInitHoles=True            # Fill holes in mask before running ACWE

# Inform user
verbose = True  # Inform user about which image is being processed

# # Time ACWE
# timeFile = os.path.join(ROOT_DIR,'MagStandardI/') + CR + '_timeStandard_FluxImb.csv'

# In[3]
# Key functions

# Updated ACWE Itteration Fuction
def itterateACWE2(I,im_size,sd_mask,m,foreground_weight=1,
                  background_weight=1/50.,inunipolarityweight=100,
                  outunipolarityweight=1/10,narrowband=2,N=10,
                  fillInitHoles=True,verbose=False):
    # Set up variables for ACWE iterations
    if fillInitHoles:
        m_seg = sp.ndimage.binary_fill_holes(m)
    else:
        m_seg = m # image to keep track of current initialization of ACWE

    I_seg = I * 1 # image of current segmentation
    if len(I.shape) == 2:
        I_seg[~sd_mask] = I[~m_seg&sd_mask].mean() # set pixels outside SD to
                                                   # mean of background to 
                                                   # force ACWE to ignore
    else:
        for i in range(I.shape[2]):
            I_segMini = I_seg[:,:,i]
            I_segMini[~sd_mask] = I[~m_seg&sd_mask].mean()

    counter = 0 # to keep track of proxy of iterations
    seg_diff_cum = np.zeros(m.shape) # to keep track of how many times
                                     # pixels change classes over
                                     # iterations
    iterate = 1 # flag to continue iterating

    # Continue iterating with N iterations of the ACWE evolution, 
    # followed by check for convergence. Running ACWE for default of N=10
    # iterations is a good trade-off between checking too often and not
    # often enough for convergence.

    while iterate:
        seg = acwe_exp2.acwe(I_seg,m_seg,N,(0,0,foreground_weight,
                                            background_weight,
                                            inunipolarityweight,
                                            outunipolarityweight),
                             narrowband,verbose) # Evolve ACWE 
                                                 # for N
                                                 # Iterations

        # update current segmentation
        I_seg = I * 1 # image of current segmentation
        if len(I.shape) == 2:
            I_seg[~sd_mask] = I[~seg&sd_mask].mean() # set pixels outside SD to
                                                     # mean of background to 
                                                     # force ACWE to ignore
        else:
            for i in range(I.shape[2]):
                I_segMini = I_seg[:,:,i]
                I_segMini[~sd_mask] = I[~seg&sd_mask].mean()

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
        percent_new_diff = float(((seg_diff_cum==1)*abs(seg_diff)).sum())/\
                           float(((seg_diff_cum>=1)*abs(seg_diff)).sum()+\
                           np.finfo(float).eps)*100 
        if verbose:
            print(str(percent_diff) + ' ' + str(percent_new_diff))
        if percent_new_diff==0 | ~(seg.sum()>0):
            iterate = 0
        if np.sum((~seg&sd_mask).astype(int)) == 0:
            iterate = 0
    return seg

# Open EUV Image
def openAIA(filename):
    
    # Extract Image and Header Data
    hdulist = fits.open(filename)
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
        I = J*1
        H = dict(h)
        
    # Prepare Display Version
    Idsp = np.clip(I,20,2500)
    Idsp = np.log10(Idsp)
    Idsp = Idsp - np.min(Idsp)
    Idsp = Idsp/np.max(Idsp)
    
    return I,Idsp,H

# Open & Align Magnetogram
def openHMI(filename,Ieuv,Heuv):
    
    # Open Magnetogram
    hdulist = fits.open(filename)
    hdulist.verify('silentfix') #necessary for successful data read
    h_mag = hdulist[1].header
    J_mag = hdulist[1].data
    hdulist.close()

    # Create Maps
    mapHMI = sunpy.map.Map((J_mag,h_mag))
    mapHMI.plot_settings['cmap']='hmimag'
    mapACWE = sunpy.map.Map(Ieuv,Heuv)

    # Reproject Magnetogram
    hmiReproject = mapHMI.reproject_to(mapACWE.wcs)
    # hmiReproject,footprint = reproject_interp(mapHMI,mapACWE.wcs,
    #                                           mapACWE.data.shape)
    I = hmiReproject.data
    # Undo Keword Renaming
    H = dict()
    for k in hmiReproject.meta.keys(): 
        H[k.upper()] = hmiReproject.meta[k] 
    
    return I,H # hmiReproject


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

#  Find location of Seed
crSeedFolder = seedFolder + CR + '/'

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
        Itmp,_,Htmp = openAIA(dataFolder+str(CR) +'/'+files[0])
        Ishape   = np.hstack([Itmp.shape,len(files)])
        I = np.zeros(Ishape); I[:,:,0] = Itmp
        H = [];               H.append(Htmp)
        
        # Open Rest
        for j in range(1,len(files)):
            # Magnetogram
            if 'hmi.' in files[j]:
                Itmp,Htmp = openHMI(dataFolder+str(CR) +'/'+files[j], I[:,:,0], H[0])
                Itmp[np.isnan(Itmp)] = np.nanmin(Itmp)
            # Remaining EUVs
            else:
                Itmp,_,Htmp = openAIA(dataFolder+str(CR) +'/'+files[j])
            # Save
            I[:,:,j] = Itmp
            H.append(Htmp)
            
        # Inform User
        if verbose:
            print('    Getting Seed')
        
        # Get seed File Name
        seedFile = seedPrefix + os.path.basename(files[0])
        for j in range(1,len(acweChoices)):
            seedFile = seedFile + '.' + keys[acweChoices[j]]
        seedFile = crSeedFolder + acweFolder + seedFile  + '.npz'
        
        # Get Seed
        _,AH,_ = as5.openSeg(seedFile)
        m      = AH['INIT_MASK']
        
        # Shrink Seed - Morphology applied to ensure small regions are retained
        #               Resized with neariest neighbor interpolation, and
        #               anti_aliasing=False.
        m = skimage.morphology.dilation(m,np.ones([int(resize_param/2),
                                                   int(resize_param/2)]))
        m = skimage.transform.resize(m.astype(float),
                                     np.asarray(m.shape)/resize_param,
                                     order=0,preserve_range=True,
                                     anti_aliasing=False)
        m = m > 0
        
        # Get Alpha
        alpha  = AH['INIT_ALPHA']
        alphar = AH['ALPHA']
        
        # Inform user
        if verbose:
            print('    Running ACWE')
            
        # # Time
        # start = time.time()
        
        # Prepare for ACWE - 1 loop to save time
        # Declare Placeholders
        if resize_param > 1:
            Ishape[:-1] = (Ishape[:-1]/resize_param).astype(int) # Shape of resized image
        Imini      = np.zeros(Ishape)    # Resized Image Array
        sun_radius = np.zeros(Ishape[0]) # Sun Radius Data
        sun_center = np.zeros(np.asarray([Ishape[0],2])) # Sun Center Data
        sd_mask    = np.zeros(Ishape).astype(int)        # Solar Disk Masks
        for j in range(len(files)):
            # HMI
            if 'hmi.' in files[j]:
                # Resize
                Imini[:,:,j],_,sun_radius[j],sun_center[j] = af6.resize_EUV(I[:,:,j],
                                                                            H[j-1],
                                                                            resize_param)
            # EUV
            else:
                # Resize
                Imini[:,:,j],_,sun_radius[j],sun_center[j] = af6.resize_EUV(I[:,:,j],
                                                                            H[j],
                                                                            resize_param)
                # Correct for Limb Brightening
                if correctLimbBrightening:
                    Imini[:,:,j] = clb.correct_limb_brightening(Imini[:,:,j],
                                                                sun_center[j],
                                                                sun_radius[j])
                # Inital Seeding
                if rollingAlpha != 0 and alpha[j] > 0:
                    sd_mask[:,:,j],_,_ = af6.inital_masks(Imini[:,:,j],
                                                          Ishape[:-1],
                                                          sun_radius[j],
                                                          sun_center[j],
                                                          alpha[j],
                                                          rollingAlpha)
                elif alpha[j] > 0:
                    sd_mask[:,:,j],_ = af6.inital_masks(Imini[:,:,j],
                                                        Ishape[:-1],
                                                        sun_radius[j],
                                                        sun_center[j],
                                                        alpha[j],
                                                        rollingAlpha)
                    
            
        # Combine Solar mask
        sd_mask = np.sum(sd_mask,axis=2).astype(bool)
        
        # RUN ACWE
        seg = itterateACWE2(Imini,Ishape,sd_mask,m,foreground_weight,
                            background_weight,inunipolarityweight,
                            outunipolarityweight, narrowband,N,fillInitHoles,
                            acweVerbose)
        
        # # Time
        # end = time.time()
        # timeTotal = end-start
        # timeTotal = str(timeTotal)
        
        # Inform User
        if verbose:
            print('    Saving Results\n\n')
        
        # Save Result
        init_mask_method = 'alpha*mean(qs)'
        as5.saveSeg(crSaveFolder + acweFile,seg,H,correctLimbBrightening,
                    resize_param,foreground_weight,background_weight,m,
                    init_mask_method,fillInitHoles,alpha,
                    alphar,narrowband,N)
        
        # # Time
        # row = acweFile + ',' + timeTotal + '\n'
        # with open(timeFile,'a+') as f:
        #     f.write(row)
        
# In[7]
# End Process

print('**Process Complete**')
