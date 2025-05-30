#-------------------------------------------------------------------------------
# Active Contours Without Edges
#
# Expeirmental Version 1: Addition of metric for estimating unipoliarity
# 
# seg = acwe(I,m,N,weights,narrowband,plot_progress)
#
# Sets up forces for implementation of the active contours without edges (ACWE)
# algorithm by Chan & Vese [1] or Chan, Sandber & Vese [2].  Calls 
# level_set_evolve for evolution of the contour.
#
# Inputs:
#     I: input image
#     m: mask (initialization for active contours)
#     N: number of iterations
#     weights: vector of weights for different energy terms
#        [mu,nu,lambda_i,lambda_o] = [surface tension (length), stiffness 
#                                     (area), inside contour, outside contour]
#     narrowband: the number of pixels over which to evaluate energies [4]
#     plot_progress: flag for whether to plot progress each iteration
#
# Output:
#    seg: mask of segmented image, corresponding to the zero level set; '1' 
#         denotes the interior of the contour (foreground) and '0' the exterior
#         (background)
# 
# References:
# [1] T. F. Chan & L. A. Vese, "Active Contours Without Edges," IEEE 
#     Transactions on Image Processing, vol. 10, pp. 266-277, 2001.
# [2] T. F. Chan, B. Y. Sandberg, & L. A. Vese, "Active Contours without Edges
#     for  Vector-Valued Images," Journal of Visual Communication and Image
#     Representation, vol. 11, no. 2, pp. 130-141, 2000.
# [3] H.-K. Zhao, T. Chan, B. Merriman, & S. Osher, "A Variational Level Set
#     Approach to Multiphase Motion," Journal of Computational Physics, vol. 
#     127, pp. 179-195, 1996.
# [4] D. Peng, B. Merriman, S. Osher, H. Zhao, & M. Kang, "A PDE-Based Fast
#     Local Level Set Method," Journal of Computational Physics," vol. 155, 
#     pp. 410-438, 1999.
# 
# Notes:
# Ideas taken from the kevin-keraudren translation of the MATLAB creaseg code.
# Requires the following python modules: scipy, numpy, matplotlib
#
# Copyright 2019 Laura Boucheron
# This file is part of ACWE.
#
# ACWE is free software: you can redistribute it and/or modify it under the 
# terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later 
# version.
#
# ACWE is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more 
# details.
#
# You should have received a copy of the GNU General Public License along with 
# ACWE.  If not, see <https://www.gnu.org/licenses/>.
#
# Updates: Sep 06, 2022: reset phi has been updated to use segmentation
#                        (seg) for all terms
#          Sep 18, 2023: Now accepts mulispeciral images (any size) following
#                        process outlined in [2].
#          Oct 13, 2023: Image Force Test - skew

from scipy.ndimage.morphology import distance_transform_edt
from numpy import log10, array, sqrt, zeros, max, mean
from matplotlib import pyplot as plt
from scipy.ndimage.filters import convolve

def sobel_gradient(I):
    # define gradient of image using sobel masks
    sx = array([[-1, 0, 1],[-2, 0, 2],[-1, 0, 1]]) # vertical sobel mask
    sy = array([[-1, -2, -1],[0, 0, 0],[1, 2, 1]]) # horizontal sobel mask
    gx = convolve(I,sx,mode='reflect')
    gy = convolve(I,sy,mode='reflect')
    g = sqrt(gx**2 + gy**2)
    return g

def level_set_evolve(F,phi,narrowband):
    # evolve level set
    phi_grad = sobel_gradient(phi) # determine gradient of phi
    delta_t = 0.49*1/max(F) # define small enough timestep per Courant-
                             # Friedrichs-Lewy (CFL) stability [REF?]
    phi[abs(phi)<=narrowband] = phi[abs(phi)<=narrowband] - \
        delta_t*F*phi_grad[abs(phi)<=narrowband] # evolve phi only in narrowband
    return phi

def acwe(I,m,N,weights,narrowband,plot_progress):
    mu = weights[0] # surface tension (length) weight
    nu = weights[1] # stiffness (area) weight
    lambda_i = weights[2] # inside contour weight
    lambda_o = weights[3] # outside contour weight
    lambda_u = weights[4] # inside unpolarity weight
    lambda_m = weights[5] # outside unpolarity weight
    epsilon  = weights[6] # Prevent divide by zero

    # Convert initial mask into signed distance function [3]
    # By convention, distance is 0 on contour, <0 outside of contour, and >0
    # inside contour [1]
    # The distance_transform_edt function returns distance to nearest 
    # *background* pixel, so need to complement mask in these operations.
    # Also need to assure that the pixels on the interior boundary of the 
    # original mask are assigned distance -0.5 pixels and those pixels on the 
    # exterior boundary assigned distance +0.5 pixels 
    phi = distance_transform_edt(~m) - distance_transform_edt(m) + m - 0.5
    
    counter = 0 # keep track of number of iterations
    iterate = 1 # flag to keep iterating
    seg = m
    while (counter<N and iterate):
        
        # Standard (Grayscale) Image
        if len(I.shape) == 2:
            
            # Standard ACWE Image Forces
            m_i = I[seg].mean() # mean of interior
            m_o = I[~seg].mean() # mean of exterior
            F_image = -lambda_i*(I[abs(phi)<=narrowband]-m_i)**2 + \
                +lambda_o*(I[abs(phi)<=narrowband]-m_o)**2 # define image force
                
            ### New Image Force - Skew Foreground ###
            # Select Region
            segInt = I[seg].flatten()
            seglen = len(segInt)
            
            # Calculate Mean, varience and non-central moment
            segm = mean(segInt)
            segvar = segInt.var()
            segm3  = mean(segInt**3)
            
            # Estimate new mean, varience, skew, an non-central moment
            mz   = (I[abs(phi)<=narrowband] + segm * seglen)/\
                   (seglen + 1)
            varz = (segvar*seglen + (I[abs(phi)<=narrowband]-mz)\
                    *(I[abs(phi)<=narrowband]-segm))/\
                   (seglen + 1)
            stdz = sqrt(varz)
            mz3  = (I[abs(phi)<=narrowband]**3 + \
                    segm3 * seglen)/\
                   (seglen + 1)
                   
            # Estimate skew
            skewz = (mz3 - 3*mz*varz - mz**3)/(stdz**3)
            
            # Apply force
            F_image = F_image -lambda_u/(abs(skewz) + epsilon)
            
            ### New Image Force - Skew Background ###
            # Select Region
            segInt = I[~seg].flatten()
            seglen = len(segInt)
            
            # Calculate Mean, varience and non-central moment
            segm = mean(segInt)
            segvar = segInt.var()
            segm3  = mean(segInt**3)
            
            # Estimate new mean, varience, standard devation, and non-central moment
            mz   = (I[abs(phi)<=narrowband] + segm * seglen)/\
                   (seglen + 1)
            varz = (segvar*seglen + (I[abs(phi)<=narrowband]-mz)\
                    *(I[abs(phi)<=narrowband]-segm))/\
                   (seglen + 1)
            stdz = sqrt(varz)
            mz3  = (I[abs(phi)<=narrowband]**3 + \
                    segm3 * seglen)/\
                   (seglen + 1)
                   
            # Estimate skew
            skewz = (mz3 - 3*mz*varz - mz**3)/(stdz**3)
            
            # Apply force
            F_image = F_image +lambda_m/(abs(skewz) + epsilon)
            
        # Vector-Valued Image [2]
        else: 
            F_image = 0
            for i in range(array(I.shape)[2]):
                
                # Standard ACWE Image Forces
                Imini = I[:,:,i]
                m_i = Imini[seg].mean()
                m_o = Imini[~seg].mean()
                F_image = -lambda_i[i]*(Imini[abs(phi)<=narrowband]-m_i)**2 + \
                          +lambda_o[i]*(Imini[abs(phi)<=narrowband]-m_o)**2 + \
                          F_image
                
                ### New Image Force - Skew Foreground ###
                # Select Region
                segInt = Imini[seg].flatten()
                seglen = len(segInt)
                
                # Calculate Mean, varience and non-central moment
                segm = mean(segInt)
                segvar = segInt.var()
                segm3  = mean(segInt**3)
                
                # Estimate new mean, varience, standard devation, and non-central moment
                mz   = (Imini[abs(phi)<=narrowband] + segm * seglen)/\
                       (seglen + 1)
                varz = (segvar*seglen + (Imini[abs(phi)<=narrowband]-mz)\
                        *(Imini[abs(phi)<=narrowband]-segm))/\
                       (seglen + 1)
                stdz = sqrt(varz)
                mz3  = (Imini[abs(phi)<=narrowband]**3 + \
                        segm3 * seglen)/\
                       (seglen + 1)
                       
                # Estimate skew
                skewz = (mz3 - 3*mz*varz - mz**3)/(stdz**3)
                
                # Apply force
                F_image = F_image -lambda_u[i]/(abs(skewz) + epsilon)
                
                ### New Image Force - Skew Background ###
                # Select Region
                segInt = Imini[~seg].flatten()
                seglen = len(segInt)
                
                # Calculate Mean, varience and non-central moment
                segm = mean(segInt)
                segvar = segInt.var()
                segm3  = mean(segInt**3)
                
                # Estimate new mean, varience, skew, an non-central moment
                mz   = (Imini[abs(phi)<=narrowband] + segm * seglen)/\
                       (seglen + 1)
                varz = (segvar*seglen + (Imini[abs(phi)<=narrowband]-mz)\
                        *(Imini[abs(phi)<=narrowband]-segm))/\
                       (seglen + 1)
                stdz = sqrt(varz)
                mz3  = (Imini[abs(phi)<=narrowband]**3 + \
                        segm3 * seglen)/\
                       (seglen + 1)
                       
                # Estimate skew
                skewz = (mz3 - 3*mz*varz - mz**3)/(stdz**3)
                
                # Apply force
                F_image = F_image +lambda_m[i]/(abs(skewz) + epsilon)
                
                
        F_length = mu*0
        F_area = nu*0

        F_total = F_length + F_area + F_image # total energy only for narrowband
    
        phi = level_set_evolve(F_total,phi,narrowband)
        seg = phi<=0

        # reset phi to signed distance transform manually 
        # could probably reinitialize phi, but code wasn't working and manual
        # approach worked
        phi = distance_transform_edt(~seg) - distance_transform_edt(seg) + \
              seg - 0.5

        if plot_progress:
            plt.figure(1)
            plt.clf()
            plt.imshow(log10(I),vmin=log10(100),vmax=log10(2500),cmap='gray')
            plt.contour(phi,0,colors='y')
            plt.axis('off')
            plt.draw()
            plt.pause(1)

        counter = counter + 1
    return seg
