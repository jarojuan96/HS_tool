#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:10:15 2026

@author: jroger
"""

import numpy as np


##### Radiance masking

def masking_asking(): #Questions for radiance thresholding
    
    mask_flag_lower_approval = False
    while mask_flag_lower_approval == False:
        mask_flag_lower = input('Masking lower values? y/n: ')
        if mask_flag_lower == 'y':
            umbral_mask_lower = float(input('Enter lower threshold: '))
            mask_flag_lower = True
            mask_flag_lower_approval = True
        elif mask_flag_lower == 'n':
            mask_flag_lower_approval = True
            mask_flag_lower = False
            umbral_mask_lower = 0
        else:
            pass
    mask_flag_upper_approval = False
    while mask_flag_upper_approval == False:
        mask_flag_upper = input('Masking upper values? y/n: ')
        if mask_flag_upper == 'y':
            umbral_mask_upper = float(input('Enter upper threshold: '))
            mask_flag_upper = True
            mask_flag_upper_approval = True
        elif mask_flag_upper == 'n':
            mask_flag_upper_approval = True
            umbral_mask_upper = 0
            mask_flag_upper = False
        else:
            pass
    return mask_flag_lower, umbral_mask_lower, mask_flag_upper, umbral_mask_upper


def mask_upper_idx(img_ref, idx, thres): #Allows for a upper radiance threshold
    
    mask = img_ref.copy()[:,:,idx] > thres
    mask = mask*1.000
    
    return mask

def mask_lower_idx(img_ref, idx, thres): #Allows for a lower radiance threshold
    
    mask = img_ref.copy()[:,:,idx] < thres
    mask = mask*1.000
    
    return mask

def mask_threshold (mask_flag_lower, thres_mask_lower, mask_flag_upper, thres_mask_upper, idx_ref, img): #provides the whole image excluding pixels according to thresholding

    if mask_flag_lower == True:
        mask_l = mask_lower_idx(img, idx_ref, thres_mask_lower)
    else:
        mask_l = np.zeros(img[:,:,idx_ref].shape)*0.000
    if mask_flag_upper == True:
        mask_u = mask_upper_idx(img, idx_ref, thres_mask_upper)
    else:
        mask_u = np.zeros(img[:,:,idx_ref].shape)*0.000

    mask_t = mask_u + mask_l
    mask_t[mask_t > 1.000] = 1.000
    img[mask_t == 1.000] = np.nan  

    return img