#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 23 12:27:31 2024

@author: jgra
"""

# In[1]
import numpy as np
from astropy.io import fits
import sunpy.map
from aiapy.calibrate import update_pointing, register

# In[2]
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

# In[3]
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