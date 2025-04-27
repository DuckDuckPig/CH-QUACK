#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 10:57:57 2023

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

    # Cycle through times
    for i in range(len(times)-1):
        
        # Time Data
        time1 = times[i]
        time2 = times[i+1]
        
        # Convert to datetime object
        timeobj1 = dmt.timeUnformat(time1)
        timeobj2 = dmt.timeUnformat(time2)
     
        # Calculate Time Difference/Gap
        gap = timeobj2 - timeobj1
        
        # Update Biggest Gap
        if i == 0:
            biggestGap = gap * 1
        elif gap > biggestGap:
            biggestGap = gap * 1
            
    print('Largest time gap in',
          os.path.basename(CR).split('.')[0],
          ':',biggestGap)
