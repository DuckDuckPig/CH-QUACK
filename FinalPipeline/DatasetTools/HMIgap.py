#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 27 13:44:22 2025

@author: jgra
"""

# In[1]
# Import libraries and tools
import pandas as pd
import DataManagmentTools as dmt
import glob
import os

# In[2]
# Dataset

traceFolder = 'DailyCadenceDownloadLists/'

# Record Time
time = 'T_REC'
mag  = 'Magnetogram'

# In[3]
# Get Dataset

CRs = sorted(glob.glob(traceFolder + '*.csv'))

# In[4]
# Cycle through dataset

for CR in CRs:

    # Open Carrington Rotation Document
    data = pd.read_csv(CR,header=0)
    keys = data.keys()
    
    # Get times
    times = data[time]
    mags  = data[mag]

    # Cycle through times
    for i in range(len(times)):
        
        # Time Data
        time1 = times[i]
        time2 = mags[i]
        
        # Convert to datetime object
        timeobj1 = dmt.timeUnformat(time1)
        timeobj2 = dmt.timeFromFilenameMag(time2)
     
        # Calculate Time Difference/Gap
        gap = timeobj2 - timeobj1
        
        # Update Biggest Gap
        if i == 0:
            biggestGap = abs(gap) * 1
        elif abs(gap) > biggestGap:
            biggestGap = abs(gap) * 1
            
    print('Largest HMI to AIA time gap in',
          os.path.basename(CR).split('.')[0],
          ':',biggestGap)
