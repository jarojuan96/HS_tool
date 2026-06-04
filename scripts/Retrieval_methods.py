#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:12:21 2026

@author: jroger
"""

import numpy as np
from scipy.signal import medfilt2d
from skimage import measure

def masking_plumes(mf, area_threshold, std):
    
    ind1 = medfilt2d(mf, kernel_size = 3) > 2*std
    ind2 = medfilt2d(mf, kernel_size = 3) < 100000  # A number that if too high, its is physically impossible.
    mask = ind1 & ind2
    
    labels = measure.label(mask, connectivity=1) 
    labels_prop = measure.regionprops_table(labels, properties=(
        'label', 'filled_area', 'eccentricity', 'major_axis_length', 'centroid', 'area'))
    labels_invalid_1 = labels_prop['label'][(labels_prop['filled_area'] < area_threshold)]
    for invalid1 in labels_invalid_1:
        labels[labels == invalid1] = 0
    
    mask = np.array(labels, dtype=bool)

    return mask


def AT_LMF(im, k_array, band_idxs_array, bool_mask, mask):

    img = im[:,:,band_idxs_array] #La imagen tiene que estar orientada/ AT==col
    k_array = k_array[:,band_idxs_array]
    M,N,B = img.shape
    LMF = np.zeros((M,N))
    img[img<=0] = np.nan

    for i in range(N):
        a, a_spar = 1, 1
        for j in range(B):
            a = a * img[:,i,j]
            if bool_mask:
                a_spar = a_spar * img[:,i,j].copy() * mask[:,i].copy()
            else:
                a_spar = a.copy()
                
        
        m1 = ~np.isnan(a)
        m2 = a != -9999.00000
        m3 = m1*m2
        idxs_notnan = np.where(m3 == True)[0]

        size_idxs = len(idxs_notnan)
        ones_array = np.ones((size_idxs,B))
        col_notnan = img[idxs_notnan,i]        
        
        idxs_notnan_spar = np.where(~np.isnan(a_spar))[0]
        col_notnan_spar = img[idxs_notnan_spar,i]
        mu = np.nanmean(np.log(col_notnan_spar),axis=0)
        mu_array = ones_array*mu
        
        alfa = 0
        
        try:
            
            cov = np.cov(np.log(col_notnan_spar), rowvar=False)
            cov_inv = np.linalg.pinv(cov)

            num = np.dot(np.dot(np.log(col_notnan) - mu_array,cov_inv),k_array[i])
            den = np.dot(np.dot(k_array[i],cov_inv),k_array[i])

            alfa = num/den
            
        except:
            pass
        
            
        LMF[idxs_notnan, i] = alfa

    return LMF

def lmf_retrieval(img, k, window, mission_name, bool_libradtran):
    
    row, col, B = img.shape
    k_arr = np.zeros((col, B))      
    k_arr[:,window] = k.copy()

    bool_mask = False; mask = None
    lmf = AT_LMF(img, k_arr, window, bool_mask, mask) 
    
    area_thres = 10 #minimum number of pixels of the positive signal cluster #Minimum detection size of plume - accoridng to Roger et al., 2025
    mask = masking_plumes(lmf, area_thres, np.nanstd(lmf))*1 #It captures the cluster meeting the condition
    mask_ = np.ones(mask.shape)*np.nan
    mask_[mask==0] = 1

    bool_mask = True; mask = mask_.copy()
    lmf = AT_LMF(img, k_arr, window, bool_mask, mask)    

    return lmf



def AT_MF(im, k_array, band_idxs_array, bool_mask, mask):

    img = im[:,:,band_idxs_array] #La imagen tiene que estar orientada/ AT==col
    k_array = k_array[:,band_idxs_array]
    M,N,B = img.shape
    MF = np.zeros((M,N))

    for i in range(N):
        a, a_spar = 1, 1
        for j in range(B):
            a = a * img[:,i,j]
            if bool_mask:
                a_spar = a_spar * img[:,i,j].copy() * mask[:,i].copy()
            else:
                a_spar = a.copy()

        m1 = ~np.isnan(a)
        m2 = a != -9999.00000
        m3 = m1*m2
        idxs_notnan = np.where(m3 == True)[0]

        size_idxs = len(idxs_notnan)
        ones_array = np.ones((size_idxs,B))
        col_notnan = img[idxs_notnan,i]        
        
        idxs_notnan_spar = np.where(~np.isnan(a_spar))[0]
        col_notnan_spar = img[idxs_notnan_spar,i]
        mu = np.nanmean(col_notnan_spar,axis=0)
        mu_array = ones_array*mu
        t = mu*k_array[i]
        
        alfa = 0
        try: #SVD error viene probablemente por una escasa estadistica (muchos NaN). Hacemos la columna = 0.
            cov = np.cov(col_notnan_spar, rowvar=False)
            cov_inv = np.linalg.pinv(cov)
        
            num = np.dot(np.dot((col_notnan-mu_array),cov_inv),t)
            den = np.dot(np.dot(t,cov_inv),t)

            alfa = num/den
            MF[idxs_notnan, i] = alfa
            
        except:
            continue

    return MF

def mf_retrieval(img, k, window, mission_name, bool_libradtran): #Here we correct for the sparsity assumption

    row, col, B = img.shape
    k_arr = np.zeros((col, B))      
    k_arr[:,window] = k.copy()

    bool_mask = False; mask = None
    mf = AT_MF(img, k_arr, window, bool_mask, mask) 
    
    area_thres = 10 #minimum number of pixels of the positive signal cluster #Minimum detection size of plume - accoridng to Roger et al., 2025
    mask = masking_plumes(mf, 10, np.nanstd(mf))*1 #It captures the cluster meeting the condition
    mask_ = np.ones(mask.shape)*np.nan
    mask_[mask==0] = 1

    bool_mask = True; mask = mask_.copy()
    mf = AT_MF(img, k_arr, window, bool_mask, mask)
    
    return mf


def AT_COMBOMF(img, k_arr_swir, window_swir, k_arr_2300, window_2300, mission_name, bool_libradtran):
    
    mf_swir = mf_retrieval(img, k_arr_swir, window_swir, mission_name, bool_libradtran)
    mf_2300 = mf_retrieval(img, k_arr_2300, window_2300, mission_name, bool_libradtran)
    
    std_swir = np.std(mf_swir)
    std_2300 = np.std(mf_2300)
    f = std_2300/std_swir

    row, col = mf_swir.shape
    
    mf_combo = mf_swir.copy()
    mf_combo_res = np.reshape(mf_combo, (row*col))
    mf_2300_res = np.reshape(mf_2300, (row*col))
    
    idxs_a = np.where(mf_combo_res > mf_2300_res)[0]
    idxs_b = np.where(mf_combo_res <= mf_2300_res)[0]
    
    #Higher score probably means new artifact - We keep 2300nm retrieval values
    mf_combo_res[idxs_a] = mf_2300_res[idxs_a]
    #We scale lower score pixels to keep plume intensity levels    
    mf_combo_res[idxs_b] = mf_combo_res[idxs_b]*f
    
    mf_combo = np.reshape(mf_combo_res, (row, col))
    
    return mf_combo, mf_2300, mf_swir
