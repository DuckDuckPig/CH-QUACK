"""
Description:
    test_acwe.py from "ACWE_python_v3" converted into callable functions.

Created on Mon Jun  8 14:28:32 2020
Edited on  Thu Apr  1 14:25:54 2021 - Adjusted for STEREO A and B
           Wed Jul  7 12:14:39 2021 - Rolling Alpha Correction,
                                      Optimized Confidence Map Implementation
           Wed Apr 27 09:59:53 2022 - Return Init Mask, Documentation
           Thu May 26 16:52:18 2022 - Return Init Mask Without Hole Filling
           Mon Feb 20 14:42:47 2023 - Condenced to tools used in paper, with 
                                      aditinal documentation
           Wed Oct 23 09:04:35 2024 - Vector Valued Functions with Unipolarity
           Fri Feb 21 13:45:49 2025 - Updated Documentation

@author: jgra
"""

# In[1]
# Import Libraries and tools
import numpy as np
from matplotlib import pyplot as plt
import skimage.transform
import copy
from .ACWE_python_v3 import correct_limb_brightening as clb
import scipy as sp
from .ACWE_python_v3 import acwe_exp2

# In[2]:
# Resizing Function
def resize_EUV(J,h,resize_param=8,interpolation='Bi-cubic'):
    '''
    Resizes solar EUV image based on user-specified resize parameter, and 
    return useful metadata about resized image.
    
    Parameters
    ----------
    J : [float]
        Solar EUV image
    h : dict
        original .fits header
    resize_param : int
        Resize Parameter which defines the degree to which the image should be
        spatially downsampled. Note that the parameter can be interpreted as 
        the denominator of a fraction such that a value of '2' would indicate 
        that the output image should be 1/2 the original spatial resolution 
        in each dimension, 4 indicates 1/4 scale, etc. By convention, this 
        parameter is 8 for SDO-AIA images and 4 for STEREO A and STEREO B 
        images, this produces an output that is 512x512 pixels in size.
    interpolation : str, optional
        Interpolation method for the resizing process. Valid options are 
        'Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
        'Bi-quartic', and 'Bi-quintic'.
        
        Default Value: 'Bi-cubic'
    Returns
    -------
        I : [float]
            Solar EUV Image, resized to user-specified dimensions.
        im_size : [int]
            Array with provides the dimensions of the image.
        sun_radius : float
            Radius of the Sun in the resized image.
        sun_center : [int]
            Coordinates of the center of the sun in the resized image.
    '''
    
    # Determine Downsample Method
    downsample = ['Nearest-neighbor','Bi-linear','Bi-quadratic','Bi-cubic',
                  'Bi-quartic','Bi-quintic']
    for order in range(len(downsample)):
        if downsample[order] == interpolation:
            break
    
    # Resize image
    if resize_param > 1:
        I = skimage.transform.resize(J,np.asarray(J.shape)/resize_param,
                                     order=order,preserve_range=True,
                                     anti_aliasing=True)
    else:
        I = copy.deepcopy(J)
        
    # Determine characteristics of image
    im_size = np.asarray(I.shape) # size of the image
    
    # Get Solar Radius from Metadata
    try: # SDO/AIA
        sun_radius = h['R_SUN']/resize_param
    except: # STERO A or B
        sun_radius = (h['RSUN']/h['CDELT1'])/resize_param
    
    # Get Solar Center from Metadata
    sun_center = np.asarray([int(round(h['CRPIX1']))-1,
                             int(round(h['CRPIX2']))-1])
    sun_center = sun_center/resize_param
    
    # Return Resized Image, Image Dimensions and Solar Radius & Center
    return I,im_size,sun_radius,sun_center


# In[3]
# Masking Functions

# Circle Mask
def make_circle_mask(c,im_dims,r):
    '''
    Defines a binary image with image dimensions im_dims of a circle with 
    center c and radius r
    
    Parameters
    ----------
    c : [int,int]
        x, and y coordinates of the center of a circle, obeying right-hand 
        rule where the upper left corner is the origin.
    im_dims : [int,int]
        Image dimensions in format [x,y] where x is the height and y is the 
        length
    r : float
        radius of circle
    Returns
    -------
    c_mask : [bool]
        Mask of size im_dims where circle of radius r, centered at c, is given
        the value of 1 and all other regions are assigned a value of 0.
    '''
    cx = c[0]
    cy = c[1]
    ix = im_dims[0]
    iy = im_dims[1]
    x,y = np.meshgrid(np.arange(-(cx),(ix-cx),1),np.arange(-(cy),(iy-cy),1))
    c_mask = (x**2+y**2)<=r**2
    return c_mask

# Initial Masks
def inital_masks(I,im_size,sun_radius,sun_center,alpha=0.3,rollingAlpha=0):
    '''
    Function returns circle mask that separates on-disk and off disk areas
    and initial mask for performing ACWE.
    
    Parameters
    ----------
    I : [float]
            Solar EUV Image, resized to user-specified dimensions.
    im_size : [int]
            Array with provides the dimensions of the image.
    sun_radius : float
        Radius of the Sun in the resized image.
    sun_center : [int]
        Coordinates of the center of the sun in the resized image.
    alpha : float, optional
        Threshold parameter as expressed in [1], alpha will be multiplied by 
        mean quiet Sun intensity (QS) to generate the threshold for the 
        initial mask. It is recommended that an alpha of 0.3 be use for AIA.
        An optimal value for STEREO is still being determined.
        
        Default Value: 0.3
    rollingAlpha : float, optional
        If no mask is produced using threshold alpha * QS, the threshold will
        be replaced with (alpha + i*rollingAlpha) * QS where 'i' is the number
        of times an empty mask was produced. This will result in a mask being 
        produced, even if no CH is present on disk. Set to 0 to disable this
        process.
        
        Default Value: 0 (Mask will always be alpha * QS)
    Returns
    -------
    sd_mask : [bool]
        Mask that separates on-disk and off-disk areas
    m : [bool]
        Initial mask, note holes are not filled prior to returning mask.
    alphar : float, optional
        The alpha parameter that was actually used to generate the mask, only 
        returned if rollingAlpha != 0
    References
    ----------
    [1] 
        L. E. Boucheron, M. Valluri, and R. T. J. McAteer, "Segmentation 
        of Coronal Holes Using Active Contours Without Edges," Solar 
        Physics, vol. 291, pp. 2353-2372, 2016.
    [2] 
        C. Verbeeck, V. Delouille, B. Mampaey, & R. De Visscher, "The 
        SPoCA-suite: Software for extraction, characterization and 
        tracking of active regions and coronal holes on EUV images," 
        Astronomy & Astrophysics, vol. 561, pp. A29, 2014.
    '''
    
    # Define solar disk mask
    sd_mask = make_circle_mask(sun_center,im_size,sun_radius)
    
    # Determine threshold value for initialization of AC as percentage of QS;  
    # estimate QS as maximum bin of histogram
    Ihist,ihist_edges = np.histogram(I[sd_mask],100) # 100 bin histogram of SD
    ihist_centers = (ihist_edges[:-1]+ihist_edges[1:])/2. # bin centers
    ihistmax = Ihist.argmax() # index of maximum bin
    initial_thresh = alpha*ihist_centers[ihistmax] # bin center of maximum bin
    m = (I<=initial_thresh)*sd_mask # initial mask
    
    if rollingAlpha != 0:
        # Save old alpha parameter for exception check
        alphar = alpha * 1
    
        # Adjust Alpha Parameter, if needed and instructed 
        # to do so to ensure valid initial mask
        while np.sum(m.astype(int)) == 0:
            alphar += 0.01 # increase alpha by 1% and regenerate initial mask
            initial_thresh = alphar*ihist_centers[ihistmax] # bin center of maximum bin
            m = (I<=initial_thresh)*sd_mask # new initial mask
        
    # Return Results
        return sd_mask,m,alphar
    else:
        return sd_mask,m

# In[4]
# Seed Prefiltering
def seed_prefilter(I,im_size,sd_mask,m,foreground_weight=1,
                  background_weight=1/50.,unipolarity_foreground_weight=0,
                  unipolarity_background_weight=0,narrowband=2,N=10,
                  minItter=50,strell=5,uniCap=0.8,fillInitHoles=True,
                  verbose=False):
    
    '''
    Runs coronal hole (CH) segmentation using active contours without edges 
    (ACWE) for a specified number of cases before removing non-unipolar regions
    
    Parameters
    ----------
    I : [float]
        Solar EUV Image, resized to user-specified dimensions, with correction
        for limb brightening, if needed.
    im_size : [int]
        Array which provides the dimensions of the image.
    sd_mask : [bool]
        Mask that separates on-disk and off-disk areas.
    m : [bool]
        Inital mask.
    foreground_weight : float or [float], optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of 
        the backround_weight compared to forground_weight.
        
        Default Value: 1
    background_weight : float or [float], optional
        Weight term for the background (quiet Sun and all remaining on-disk
        features) homogeneity within the energy functional.
        
        Default Value: 1/50.0
    unipolarity_foreground_weight : float or [float], optional
        Weight term for the foreground (CH) unipolarity within the energy
        functional. We note that for seed prefiltering ALL 
        unipolarity_forground_weight values should, ideally, be 0,
        but provide the user the option should they wish to try other values.
        
        Defalut Value: 0
    unipolarity_background_weight : float or [float], optional
        weight term for the bacground (quiet Sun and all remaining on-disk
        features) unipolarity within the energy functional. We note that for 
        seed prefiltering ALL unipolarity_background_weight values should, 
        ideally, be 0,but provide the user the option should they wish to try 
        other values.
        
        Default Value : 0
    narrowband : int, optional
        Constraint on ACWE evolution to ensure iterative optimization process
        does not result in overcorrection of contour boundary.
        
        Default Value: 2
    N : int, optional
        Number of iterations of ACWE between checks for convergence
        
        Default Value: 10
    minItter : int, optional
        Minimum number of itterations to perform before prefiltering the seed.
        We note that the true number of itterations will be 
        N * np.ceil(minItter/N) times.
        
        Default Value: 50
    strell : int, optional
        Size of structuring element used to group regions when pre-filtering
        intial seed. For consistanty, this strell should be 5 for a 512x512 
        pixel image, and 10 for a 1024x1024 pixel image, and 40 for a 4096x4096
        pixel image.
        
        Default Value: 5
    uniCap : float, optional
        Maximum line-of-sight unipolarity of region, using equation described
        in [2], in a seed region, after eveovling for N * ceil(minItter/N)
        itterations. All seeds with a value > unicap will be removed.
        
        Default Value: 0.80
    fillInitHoles : bool, optional
        Fill holes in initial mask using morphology prior to performing ACWE
        
        Default Value: True
    verbose : bool, optional
        Display ACWE evolution in real time.
        
        Default Value: False
    Returns
    -------
    m_seg : [bool]
        prefiltered seed, in dimensions 
        np.asarray(J.shape)/resize_param
    References
    ----------
    [1] 
        L. E. Boucheron, M. Valluri, and R. T. J. McAteer, "Segmentation 
        of Coronal Holes Using Active Contours Without Edges," Solar
        Physics, vol. 291, pp. 2353-2372, 2016.
    [2]
        Y.-K. Ko, K. Muglach, Y.-M. Wang, P. R. Young, & S. T. Lepri, "Temporal
        evolution of solar wind ion composition and thier source coronal holes 
        during the declining phase of cycle 23. i. low-latitude extension of
        polar coronal holes," The Astrophysical Journal, vol. 787, no. 2, 
        p. 121, 2014.
    '''
    
    #### Function Prep ####
    
    # Set up variables for ACWE iterations
    if fillInitHoles:
        m_seg = sp.ndimage.binary_fill_holes(m)
    else:
        m_seg = m # image to keep track of current initialization of ACWE
    
    # Mask Off-Disk Regions
    I_seg = I * 1 # image of current segmentation
    if len(I.shape) == 2: # Single Image
        I_seg[~sd_mask] = I[~m_seg&sd_mask].mean() # set pixels outside SD to
                                                   # mean of background to 
                                                   # force ACWE to ignore
    else: # Vector-Valued Image
        for i in range(I.shape[2]):
            I_segMini = I_seg[:,:,i]
            I_segMini[~sd_mask] = I[~m_seg&sd_mask].mean()
            
    # Prepare for Itteration
    counter = 0 # to keep track of proxy of iterations
    itterCounter = int(np.ceil(minItter/N))
    
    
    # Itterate for N * ceil(switch/N) using homogenity alone, before
    # removing regions of low unipolarity
    while counter < itterCounter:
        
        seg = acwe_exp2.acwe(I_seg,m_seg,N,(0,0,foreground_weight,
                                            background_weight,
                                            unipolarity_foreground_weight,
                                            unipolarity_background_weight),
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
        
        m_seg = seg # update m_seg image
        counter = counter + 1 # iterate counter
    
    #### Remove Low-Unipolairty Regions ####
    
    # Find and Locate Segmentation Clusters
    SegClusters = seg * 1
    SegClusters = skimage.morphology.dilation(SegClusters.astype(bool),
                                              np.ones([strell,strell]))
    SegClusters = skimage.measure.label(SegClusters,connectivity=2)
    
    # Check Unipolarity of Each Region
    for i in range(1,np.max(SegClusters)+1):
        
        # Define Cluster
        cluster = SegClusters == i
        CH_kept = cluster * seg
        
        hmiReproject = I[:,:,-1]
        
        # Flatten Selected Region
        flattened = hmiReproject[CH_kept].flatten()
       
        # Remove NAN Values, if any Persist
        flattened = flattened[np.logical_not(np.isnan(flattened))]
    
        # Calculate Unweighted Unipolarity
        num = np.mean(np.abs(flattened)) - np.abs(np.mean(flattened))
        den = np.mean(np.abs(flattened))
        uni = num/den
        
        # Remove low-unipolarity regions
        if uni > uniCap:
            m_seg[CH_kept] = 0
    
    # Return pre-filtered Seed
    return m_seg

# In[4]
# ACWE Iterations

# Single ACWE Segmentation
def itterate_acwe(I,im_size,sd_mask,m,foreground_weight=1,
                  background_weight=1/50.,unipolarity_foreground_weight=1,
                  unipolarity_background_weight=1/50.,narrowband=2,N=10,
                  fillInitHoles=True,verbose=False):
    '''
    Runs coronal hole (CH) segmentation using active contours without edges 
    (ACWE) as described in [1]. This version is capable of Quantifying 
    Unipolarity via Active Contour Kinetics (QUACK), generating a single 
    segmentation based on the two goals of maximizing unipolarity and 
    maximizing homogeneity of both the foreground and background as expressed 
    in one or more channels of a vector-valued image.
    
    Parameters
    ----------
    I : [float]
        Solar EUV Image, resized to user-specified dimensions, with correction
        for limb brightening, if needed.
    im_size : [int]
        Array which provides the dimensions of the image.
    sd_mask : [bool]
        Mask that separates on-disk and off-disk areas.
    m : [bool]
        Inital mask.
    foreground_weight : float or [float], optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of 
        the backround_weight compared to forground_weight.
        
        Default Value: 1
    background_weight : float or [float], optional
        Weight term for the background (quiet Sun and all remaining on-disk
        features) homogeneity within the energy functional.
        
        Default Value: 1/50.0
    unipolarity_foreground_weight : float or [float], optional
        Weight term for the foreground (CH) unipolarity within the energy
        functional.
        
        Defalut Value: 1
    unipolarity_background_weight : float or [float], optional
        weight term for the bacground (quiet Sun and all remaining on-disk
        features) unipolarity within the energy functional.
        
        Default Value : 1.50.0
    narrowband : int, optional
        Constraint on ACWE evolution to ensure iterative optimization process
        does not result in overcorrection of contour boundary.
        
        Default Value: 2
    N : int, optional
        Number of iterations of ACWE between checks for convergence
        
        Default Value: 10
    fillInitHoles : bool, optional
        Fill holes in initial mask using morphology prior to performing ACWE
        
        Default Value: True
    verbose : bool, optional
        Display ACWE evolution in real time.
        
        Default Value: False
    Returns
    -------
    seg : [bool]
        final segmentation mask, in dimensions 
        np.asarray(J.shape)/resize_param
    References
    ----------
    [1] 
        L. E. Boucheron, M. Valluri, and R. T. J. McAteer, "Segmentation 
        of Coronal Holes Using Active Contours Without Edges," Solar
        Physics, vol. 291, pp. 2353-2372, 2016.
    [2]
        Y.-K. Ko, K. Muglach, Y.-M. Wang, P. R. Young, & S. T. Lepri, "Temporal
        evolution of solar wind ion composition and thier source coronal holes 
        during the declining phase of cycle 23. i. low-latitude extension of
        polar coronal holes," The Astrophysical Journal, vol. 787, no. 2, 
        p. 121, 2014.
    '''
    
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
                                            unipolarity_foreground_weight,
                                            unipolarity_background_weight),
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
        if np.sum(seg.astype(int)) == 0:
            iterate = 0
    return seg

# In[5]
# Running ACWE
def run_acwe(J,h,files,resize_param=8,foreground_weight=1,background_weight=1/50.,
             unipolarity_foreground_weight=1,unipolarity_background_weight=1/50.,
             alpha=0.3,narrowband=2,N=10,prefilter_foreground_weight=1,
             prefilter_background_weight=1/50.,verbose=False,
             prefilter_unipolarity_foreground_weight = 0,
             prefilter_unipolarity_background_weight = 0, minItter=50, strell=5,
             uniCap=0.8, prefilter_fillInitHoles=True, correctLimbBrightening=True,
             rollingAlpha=0,fillInitHoles=False,prefilterSeed=True):
    
    '''
    Primary function for running coronal hole (CH) segmentation using active 
    contours without edges (ACWE). This version is capable of Quantifying 
    Unipolarity via Active Contour Kinetics (QUACK), generating a single 
    segmentation based on the two goals of maximizing unipolarity and 
    maximizing homogeneity of both the foreground and background as expressed 
    in one or more channels of a vector-valued image.
    
    Parameters
    ----------
    J : [float]
        Solar EUV image stored as a numpy array
    h : dict
        .fits header for Solar EUV image J
    resize_param : int, optional
        The factor by which the image will be downsampled. In general
        operation ACWE usually operates over a 512X512 pixel image,
        thus it is recommended that 
        resize_param = np.min(np.asarray(J.shape)/np.array([512,512])).astype(int)
        
        Default Value: 8
    foreground_weight : float or [float], optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of 
        the backround_weight compared to forground_weight.
        
        Default Value: 1
    background_weight : float or [float], optional
        Weight term for the background (quiet Sun and all remaining on-disk
        features) homogeneity within the energy functional.
        
        Default Value: 1/50.0
    unipolarity_foreground_weight : float or [float], optional
        Weight term for the foreground (CH) unipolarity within the energy
        functional.
        
        Default Value: 1
    unipolarity_background_weight : float or [float], optional
        weight term for the bacground (quiet Sun and all remaining on-disk
        features) unipolarity within the energy functional.
        
        Default Value: 1/50.
    alpha : float, optional
        Threshold parameter alpha will be multiplied by mean quiet sun 
        intensity to generate the threshold for the initial mask.

        Default Value: 0.3
    narrowband : int, optional
        Constraint on ACWE evolution to ensure iterative optimization process 
        does not result in overcorrection of contour boundary.
        
        Default Value: 2
    N : int, optional
        Number of iterations of ACWE between checks for convergence.
        
        Default Value: 10
    prefilter_foreground_weight : float or [float], optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of 
        the backround_weight compared to forground_weight.
        
        Default Value: 1
    prefilter_background_weight : float or [float], optional
        Weight term for the background (quiet Sun and all remaining on-disk
        features) homogeneity within the energy functional.
        
        Default Value: 1/50.0
    verbose : bool, optional
        Display ACWE evolution in real time.
        
        Default Value: False
    prefilter_unipolarity_foreground_weight : float or [float], optional
        Weight term for the foreground (CH) unipolarity within the energy
        functional. We note that for seed prefiltering ALL 
        unipolarity_forground_weight values should, ideally, be 0,
        but provide the user the option should they wish to try other values.
        
        Defalut Value: 0
    prefilter_unipolarity_background_weight : float or [float], optional
        weight term for the bacground (quiet Sun and all remaining on-disk
        features) unipolarity within the energy functional. We note that for 
        seed prefiltering ALL unipolarity_background_weight values should, 
        ideally, be 0,but provide the user the option should they wish to try 
        other values.
        
        Default Value : 0
    minItter : int, optional
        Minimum number of itterations to perform before prefiltering the seed.
        We note that the true number of itterations will be 
        N * np.ceil(minItter/N) times.
        
        Default Value: 50
    strell : int, optional
        Size of structuring element used to group regions when pre-filtering
        intial seed. For consistanty, this strell should be 5 for a 512x512 
        pixel image, and 10 for a 1024x1024 pixel image, and 40 for a 4096x4096
        pixel image.
        
        Default Value: 5
    uniCap : float, optional
        Maximum line-of-sight unipolarity of region, using equation described
        in [2], in a seed region, after eveovling for N * ceil(minItter/N)
        itterations. All seeds with a value > unicap will be removed.
        
        Default Value: 0.80
    correctLimbBrightening : bool, optional
        Perform limb brightening correction of [2]. Note, process is not 
        recommended on STEREO images.
        
        Default Value: True
    rollingAlpha : float, optional
        If no mask is produced using threshold alpha * QS, the threshold will
        be replaced with (alpha + i*rollingAlpha) * QS where 'i' is the number
        of times an empty mask was produced. This will result in a mask being 
        produced, even if no CH is present on disk. Set to 0 to disable this
        process.
        
        Default Value: 0 (Mask will always be alpha * QS)
    fillInitHoles : bool, optional
        Fill holes in initial mask
        
        Default Value: True
        
    Returns
    -------
    seg : [bool]
        final segmentation mask, in dimensions 
        np.asarray(J.shape)/resize_param
    alphar : float, optional
        the final alpha parameter, returned if (and only if)
        oldThreshold == True and rollingAlpha == True
    m : [bool]
        Initial mask without any holes filled
    
    References
    ----------
    [1] 
        L. E. Boucheron, M. Valluri, and R. T. J. McAteer, "Segmentation 
        of Coronal Holes Using Active Contours Without Edges," Solar 
        Physics, vol. 291, pp. 2353-2372, 2016.
    [2] 
        C. Verbeeck, V. Delouille, B. Mampaey, & R. De Visscher, "The 
        SPoCA-suite: Software for extraction, characterization and 
        tracking of active regions and coronal holes on EUV images," 
        Astronomy & Astrophysics, vol. 561, pp. A29, 2014.
    '''
    Jshape = np.array(J.shape)
    
    if len(Jshape) == 2:
        # Resize image
        I,Jshape,sun_radius,sun_center = resize_EUV(J,h,resize_param)
        
        # Correct limb brightening per Verbeeck et al. 2014
        if correctLimbBrightening:
            I = clb.correct_limb_brightening(I,sun_center,sun_radius)
    
        #  Define solar disk mask and initial mask
        if rollingAlpha != 0:
            sd_mask,m,alphar = inital_masks(I,Jshape,sun_radius,sun_center,
                                            alpha,rollingAlpha)
        else:
            sd_mask,m = inital_masks(I,Jshape,sun_radius,sun_center,alpha,
                                     rollingAlpha)
    else:
        # Prepare for ACWE - 1 loop to save time
        # Declare Placeholders
        if resize_param > 1:
            Jshape[:-1] = (Jshape[:-1]/resize_param).astype(int) # Shape of resized image
        I          = np.zeros(Jshape)    # Resized Image Array
        sun_radius = np.zeros(Jshape[0]) # Sun Radius Data
        sun_center = np.zeros(np.asarray([Jshape[0],2])) # Sun Center Data
        sd_mask    = np.zeros(Jshape).astype(int)        # Solar Disk Masks
        m          = np.zeros(Jshape).astype(int)        # Inital Masks
        alphar     = np.zeros(len(alpha))   # Alpha out
        for j in range(len(files)):
            # HMI
            if 'hmi.' in files[j]:
                # Resize
                I[:,:,j],_,sun_radius[j],sun_center[j] = resize_EUV(J[:,:,j],
                                                                    h[j-1],
                                                                    resize_param)
                # Mark that not used for seeding
                alphar[j] = -1
            # EUV
            else:
                # Resize
                I[:,:,j],_,sun_radius[j],sun_center[j] = resize_EUV(J[:,:,j],
                                                                    h[j],
                                                                    resize_param)
                # Correct for Limb Brightening
                if correctLimbBrightening:
                    I[:,:,j] = clb.correct_limb_brightening(I[:,:,j],
                                                            sun_center[j],
                                                            sun_radius[j])
                # Inital Seeding
                if rollingAlpha != 0 and alpha[j] > 0:
                    sd_mask[:,:,j],m[:,:,j],alphar[j] = inital_masks(I[:,:,j],
                                                                     Jshape[:-1],
                                                                     sun_radius[j],
                                                                     sun_center[j],
                                                                     alpha[j],
                                                                     rollingAlpha)
                elif alpha[j] > 0:
                    sd_mask[:,:,j],m[:,:,j] = inital_masks(I[:,:,j],
                                                           Jshape[:-1],
                                                           sun_radius[j],
                                                           sun_center[j],
                                                           alpha[j],
                                                           rollingAlpha)
                    alphar[j] = alpha[j] * 1
            
        # Combine masks and Seeds
        sd_mask = np.sum(sd_mask,axis=2).astype(bool)
        m       = np.sum(m,axis=2).astype(bool)
        
    # Prefilter Seeed
    if prefilterSeed:
        m2 = seed_prefilter(I,Jshape,sd_mask,m,prefilter_foreground_weight,
                          prefilter_background_weight,
                          prefilter_unipolarity_foreground_weight,
                          prefilter_unipolarity_background_weight,
                          narrowband,N,minItter,strell,uniCap,
                          prefilter_fillInitHoles,verbose)
        m = np.stack([m,m2],axis=0)
    else:
        m2 = m * 1
    
    # RUN ACWE
    if np.sum(m2.astype(int)) != 0:
        seg = itterate_acwe(I,Jshape,sd_mask,m2,foreground_weight,
                            background_weight,unipolarity_foreground_weight,
                            unipolarity_background_weight,narrowband,N,
                            fillInitHoles,verbose)
    else:
        seg = m2 * 1

    
    
    # Return Results
    if rollingAlpha != 0:
        return seg,alphar,m
    
    else:
        return seg,m

# ACWE Confidence Map
def run_acwe_confidence_map(J,h,files,resize_param=8,foreground_weights=[1],
                            background_weights=[1/50.],
                            unipolarity_foreground_weights=[1],
                            unipolarity_background_weights=[1/50.],alpha=0.3,
                            narrowband=2,N=10,prefilter_foreground_weight=1,
                            prefilter_background_weight=1/50.,verbose=False,
                            prefilter_unipolarity_foreground_weight = 0,
                            prefilter_unipolarity_background_weight = 0, 
                            minItter=50, strell=5,uniCap=0.8, 
                            prefilter_fillInitHoles=True, 
                            correctLimbBrightening=True,rollingAlpha=0,
                            fillInitHoles=False,prefilterSeed=True):
    
    '''
    Primary function for running coronal hole (CH) segmentation using active 
    contours without edges (ACWE) in order to generate confidence maps. This 
    version is capable of Quantifying Unipolarity via Active Contour 
    Kinetics (QUACK), generating segmentations based on the two goals of 
    maximizing unipolarity and maximizing homogeneity of both the foreground 
    and background as expressed in one or more channels of a vector-valued 
    image.
    
    Parameters
    ----------
    J : [float]
        Solar EUV image stored as a numpy array
    h : dict
        .fits header for Solar EUV image J
    resize_param : int, optional
        The factor by which the image will be downsampled. In general
        operation ACWE usually operates over a 512X512 pixel image,
        thus it is recommended that 
        resize_param = np.min(np.asarray(J.shape)/np.array([512,512])).astype(int)
        
        Default Value: 8
    foreground_weights : [float], optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of 
        the backround_weight compared to forground_weight.
        
        Default Value: [1]
    background_weights : [float], optional
        Weight term for the background (quiet Sun and all remaining on-disk
        features) homogeneity within the energy functional.
        
        Default Value: [1/50.0]
    unipolarity_foreground_weights : [float], optional
        Weight term for the foreground (CH) unipolarity within the energy
        functional.
        
        Default Value: [1]
    unipolarity_background_weights : [float], optional
        weight term for the bacground (quiet Sun and all remaining on-disk
        features) unipolarity within the energy functional.
        
        Default Value: [1/50.0]
    alpha : float, optional
        Threshold parameter alpha will be multiplied by mean quiet sun 
        intensity to generate the threshold for the initial mask.

        Default Value: 0.3
    narrowband : int, optional
        Constraint on ACWE evolution to ensure iterative optimization process 
        does not result in overcorrection of contour boundary.
        
        Default Value: 2
    N : int, optional
        Number of iterations of ACWE between checks for convergence.
        
        Default Value: 10
    prefilter_foreground_weight : float or [float], optional
        Weight term for the foreground (CH) homogeneity within the energy 
        functional. It is recommended that this value be kept at 1 to
        facilitate an intuitive understanding of the relative strength of 
        the backround_weight compared to forground_weight.
        
        Default Value: 1
    prefilter_background_weight : float or [float], optional
        Weight term for the background (quiet Sun and all remaining on-disk
        features) homogeneity within the energy functional.
        
        Default Value: 1/50.0
    verbose : bool, optional
        Display ACWE evolution in real time.
        
        Default Value: False
    prefilter_unipolarity_foreground_weight : float or [float], optional
        Weight term for the foreground (CH) unipolarity within the energy
        functional. We note that for seed prefiltering ALL 
        unipolarity_forground_weight values should, ideally, be 0,
        but provide the user the option should they wish to try other values.
        
        Defalut Value: 0
    prefilter_unipolarity_background_weight : float or [float], optional
        weight term for the bacground (quiet Sun and all remaining on-disk
        features) unipolarity within the energy functional. We note that for 
        seed prefiltering ALL unipolarity_background_weight values should, 
        ideally, be 0,but provide the user the option should they wish to try 
        other values.
        
        Default Value : 0
    minItter : int, optional
        Minimum number of itterations to perform before prefiltering the seed.
        We note that the true number of itterations will be 
        N * np.ceil(minItter/N) times.
        
        Default Value: 50
    strell : int, optional
        Size of structuring element used to group regions when pre-filtering
        intial seed. For consistanty, this strell should be 5 for a 512x512 
        pixel image, and 10 for a 1024x1024 pixel image, and 40 for a 4096x4096
        pixel image.
        
        Default Value: 5
    uniCap : float, optional
        Maximum line-of-sight unipolarity of region, using equation described
        in [2], in a seed region, after eveovling for N * ceil(minItter/N)
        itterations. All seeds with a value > unicap will be removed.
        
        Default Value: 0.80
    correctLimbBrightening : bool, optional
        Perform limb brightening correction of [2]. Note, process is not 
        recommended on STEREO images.
        
        Default Value: True
    rollingAlpha : float, optional
        If no mask is produced using threshold alpha * QS, the threshold will
        be replaced with (alpha + i*rollingAlpha) * QS where 'i' is the number
        of times an empty mask was produced. This will result in a mask being 
        produced, even if no CH is present on disk. Set to 0 to disable this
        process.
        
        Default Value: 0 (Mask will always be alpha * QS)
    fillInitHoles : bool, optional
        Fill holes in initial mask
        
        Default Value: True
        
    Returns
    -------
    Segs : [bool]
        final segmentation mask, in dimensions 
        np.hstack([len(background_weight),np.asarray(J.shape)/resize_param])
    alphar : float, optional
        the final alpha parameter, returned if (and only if)
        oldThreshold == True and rollingAlpha == True
    m : [bool]
        Initial mask without any holes filled
    
    References
    ----------
    [1] 
        L. E. Boucheron, M. Valluri, and R. T. J. McAteer, "Segmentation 
        of Coronal Holes Using Active Contours Without Edges," Solar 
        Physics, vol. 291, pp. 2353-2372, 2016.
    [2] 
        C. Verbeeck, V. Delouille, B. Mampaey, & R. De Visscher, "The 
        SPoCA-suite: Software for extraction, characterization and 
        tracking of active regions and coronal holes on EUV images," 
        Astronomy & Astrophysics, vol. 561, pp. A29, 2014.
    '''
    Jshape = np.array(J.shape)
    
    if len(Jshape) == 2:
        # Resize image
        I,Jshape,sun_radius,sun_center = resize_EUV(J,h,resize_param)
        
        # Correct limb brightening per Verbeeck et al. 2014
        if correctLimbBrightening:
            I = clb.correct_limb_brightening.correct_limb_brightening(I,sun_center,
                                                                      sun_radius)
    
        #  Define solar disk mask and initial mask
        if rollingAlpha != 0:
            sd_mask,m,alphar = inital_masks(I,Jshape,sun_radius,sun_center,
                                            alpha,rollingAlpha)
        else:
            sd_mask,m = inital_masks(I,Jshape,sun_radius,sun_center,alpha,
                                     rollingAlpha)
    else:
        # Prepare for ACWE - 1 loop to save time
        # Declare Placeholders
        if resize_param > 1:
            Jshape[:-1] = (Jshape[:-1]/resize_param).astype(int) # Shape of resized image
        I          = np.zeros(Jshape)    # Resized Image Array
        sun_radius = np.zeros(Jshape[0]) # Sun Radius Data
        sun_center = np.zeros(np.asarray([Jshape[0],2])) # Sun Center Data
        sd_mask    = np.zeros(Jshape).astype(int)        # Solar Disk Masks
        m          = np.zeros(Jshape).astype(int)        # Inital Masks
        alphar     = np.zeros(len(alpha))   # Alpha out
        for j in range(len(files)):
            # HMI
            if 'hmi.' in files[j]:
                # Resize
                I[:,:,j],_,sun_radius[j],sun_center[j] = resize_EUV(J[:,:,j],
                                                                    h[j-1],
                                                                    resize_param)
                # Mark that not used for seeding
                alphar[j] = -1
            # EUV
            else:
                # Resize
                I[:,:,j],_,sun_radius[j],sun_center[j] = resize_EUV(J[:,:,j],
                                                                    h[j],
                                                                    resize_param)
                # Correct for Limb Brightening
                if correctLimbBrightening:
                    I[:,:,j] = clb.correct_limb_brightening(I[:,:,j],
                                                            sun_center[j],
                                                            sun_radius[j])
                # Inital Seeding
                if rollingAlpha != 0 and alpha[j] > 0:
                    sd_mask[:,:,j],m[:,:,j],alphar[j] = inital_masks(I[:,:,j],
                                                                     Jshape[:-1],
                                                                     sun_radius[j],
                                                                     sun_center[j],
                                                                     alpha[j],
                                                                     rollingAlpha)
                elif alpha[j] > 0:
                    sd_mask[:,:,j],m[:,:,j] = inital_masks(I[:,:,j],
                                                           Jshape[:-1],
                                                           sun_radius[j],
                                                           sun_center[j],
                                                           alpha[j],
                                                           rollingAlpha)
                    alphar[j] = alpha[j] * 1
            
        # Combine masks and Seeds
        sd_mask = np.sum(sd_mask,axis=2).astype(bool)
        m       = np.sum(m,axis=2).astype(bool)
        
    # Prefilter Seeed
    if prefilterSeed:
        m2 = seed_prefilter(I,Jshape,sd_mask,m,prefilter_foreground_weight,
                          prefilter_background_weight,
                          prefilter_unipolarity_foreground_weight,
                          prefilter_unipolarity_background_weight,
                          narrowband,N,minItter,strell,uniCap,
                          prefilter_fillInitHoles,verbose)
        m = np.stack([m,m2],axis=0)
    else:
        m2 = m*1        
    
    # RUN ACWE - Confidence Map
    outputShape = np.asarray(m2.shape)
    outputShape = np.hstack([len(background_weights),outputShape]).astype(int)
    Segs = np.empty(outputShape); Segs[:] = np.nan
    
    if np.sum(m2.astype(int)) != 0:
        for i in range(len(background_weights)):
            seg = itterate_acwe(I,Jshape,sd_mask,m2,foreground_weights[i],
                                background_weights[i],
                                unipolarity_foreground_weights[i],
                                unipolarity_background_weights[i],narrowband,N,
                                fillInitHoles,verbose)
            Segs[i] = seg
    else:
        Segs[:] = 0
    
    
    # Return Results
    if rollingAlpha != 0:
        return Segs,alphar,m
    
    else:
        return Segs,m
