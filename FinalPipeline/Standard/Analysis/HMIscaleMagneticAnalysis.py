#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May 11 16:20:04 2025

@author: jgra
"""

# In[1]
# Import Libaries and Tools
import os
import sys
import glob
from astropy.io import fits
import sunpy.map
import numpy as np
import skimage.morphology
import skimage.measure
import astropy.units as u

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
saveName = '.Mag2.npz'

strell = 40

# Inform User
verbose = True

# Overwrite existing data
overwrite = False

# In[3]
# List of data

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
        sMmean  = []
        usMmean = []
        cMass   = []
        
        # Cycyle Through CR
        for k in range(len(files)):
            
            file = files[k]
               
            # Inform User
            if verbose:
                print('   ',k+1,'of',len(files),'-', os.path.basename(file))
                print('        Opening File')
            
            # Open Segmentation
            H,AH,SEG = as5.openSeg(file)
            
            # Inform User
            if verbose:
                print('        Restoring Scale')
            
            # Resize File
            SEG = ars.upscale(SEG,AH)
            
            # Inform User
            if verbose:
                print('        Getting Magnetogram Data')
            
            # Find Magnetogram 
            folder  = CR + '/' + file.split('/')[-2] + '/'
            hmiFile = glob.glob(dataFolder + folder + 'hmi.*.fits')[0]
            
            # Open Magentogram
            hdulist = fits.open(hmiFile)
            hdulist.verify('silentfix') #necessary for successful data read
            Hhmi = hdulist[1].header
            Ihmi = hdulist[1].data
            hdulist.close()
            
            # Inform User
            if verbose:
                print('        Gathering Coordinate Data')
            
            # Create Map
            mapHMI = sunpy.map.Map((Ihmi,Hhmi))
            mapHMI.plot_settings['cmap']='hmimag'
                
            # Get Coordinate Information
            cord = sunpy.map.all_coordinates_from_map(mapHMI)
            
            # Convert to heliocentric
            cord = cord.transform_to('heliocentric')
            
            # Corridante Transform (Spherical)
            theta = np.arctan((np.sqrt(cord.x**2+cord.y**2)/cord.z))
            
            # Inform User
            if verbose:
                print('        Calculating baseline Area and Flux')
            
            # Area in Mm
            area = 0.132/np.cos(theta)
            
            # Weighted Flux
            flux = Ihmi/np.cos(theta)
            
            # Inform User
            if verbose:
                print('        Identifying Regions')
            
            # Find and Label CH Regions
            ACWE_clusters = SEG * 1
            ACWE_clusters[np.where(ACWE_clusters>0)] = 1
            ACWE_clusters = skimage.morphology.dilation(ACWE_clusters.astype(bool),
                                                        np.ones([strell,strell]))
            ACWE_clusters = skimage.measure.label(ACWE_clusters,connectivity=2)
            
            # Inform User
            if verbose:
                print('        Calculating Regions:')
            
            # Split CH Regions
            for j in range(1,np.max(ACWE_clusters)+1):
                
                # Inform User
                if verbose:
                    print('           ',j,'of',np.max(ACWE_clusters))
                    print('                Getting Region Data')
                
                # Convert cluster back to region
                cluster = ACWE_clusters * 1
                cluster[np.where(cluster!=j)] = 0
                cluster = cluster/np.max(cluster)
                SEG2 = SEG*cluster
                
                # Inform User
                if verbose:
                    print('                Aligning to HMI Observation')
                
                # align to HMI Observation
                mapQUACK = sunpy.map.Map(SEG2,H[0])
                mapCster = sunpy.map.Map(cluster,H[0])
                quackRep = mapQUACK.reproject_to(mapHMI.wcs)
                csterRep = mapCster.reproject_to(mapHMI.wcs)
                nSEG2    = quackRep.data; nSEG2[np.isnan(nSEG2)] = 0
                nCluster = (csterRep.data > 0).astype(int)
                
                # Inform User
                if verbose:
                    print('                Calculating Core Data')
                
                # Date
                dates.append(dmt.timeUnformat(H[0]['T_REC']))
                
                # Areas
                pxAreas.append(np.nansum(nSEG2.astype(int)))
                Ai = area*nSEG2.astype(int)
                MmAreas.append(np.nansum(Ai))
                
                # Magnetic Field
                Bi    = flux*nSEG2.astype(int)
                Phi   = np.nansum(Bi*Ai)
                Phi_t = np.nansum(np.abs(Bi)*Ai)
                sMmean.append(Phi/np.nansum(Ai))
                usMmean.append(Phi_t/np.nansum(Ai))
                
                # Get Location Data
                ACWE_cluser_Props = skimage.measure.regionprops(nCluster,nSEG2)
                centroid = ACWE_cluser_Props[0].centroid_weighted
                cMass.append(mapHMI.wcs.pixel_to_world(centroid[0]*u.pixel,
                                                       centroid[1]*u.pixel))
                
        # Save Data
        np.savez_compressed(results, dates, pxAreas, MmAreas, 
                            sMmean, usMmean, cMass)

# In[5]
# End Process

print('**Process Complete**')