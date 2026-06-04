# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import os
import numpy as np
from scipy.signal import medfilt2d
import matplotlib.pyplot as plt

from EnMAP_reader import read_enmap_img
from EMIT_reader import read_EMIT
from PRISMA_reader import read_prisma_img
from GF_reader import read_gf_img

from envi_files import simple_bsq_save_md, load_hdr_file
from radiance_masking import mask_threshold, masking_asking
from Retrieval_methods import mf_retrieval
from RT_functions import read_luts_libradtran, read_luts_libradtran_h2o, get_k_libradtran, convolution_khr

import variable_definition

import warnings
warnings.filterwarnings("ignore")



### Radiance preview for radiance thresholding

def rad_preview(p, n, bool_mask, mask_flag_lower, thres_mask_lower, mask_flag_upper, thres_mask_upper, mission_name):

    if mission_name == 'EnMAP':
        swir_flag = True #Our spectral range of interest
        img, wvl_simple, fwhm_simple, hsf, sza, vza = read_enmap_img(p, n, swir_flag)

    elif mission_name == 'EMIT':
        img, wvl_simple, fwhm_simple, sza, vza = read_EMIT(p, n)

    elif mission_name == 'PRISMA':
        swir_flag = True #Our spectral range of interest
        img, wvl_mat, fwhm_mat, vza, sza = read_prisma_img(p, n, swir_flag)
        wvl_simple = wvl_mat[500] #the central column (out of 1000)
        
    elif mission_name == 'GF5':
        img, wvl_mat, fwhm_mat, vza, sza = read_gf_img(p, n)
        wvl_simple = wvl_mat[500] 
        

    wvl_ref = 2100 #nm (not influenced by CH4, but similar radiance levels)
    idx_ref = np.abs(wvl_simple - wvl_ref).argmin()

    cmap = plt.cm.gray
    cmap.set_bad(color='darkolivegreen')    
    if bool_mask == False:
        fig = plt.figure(figsize=(10,10)); plt.imshow(img[:,:,idx_ref], cmap=cmap); ax = plt.gca(); ax.axis('off'); cbar = plt.colorbar(orientation='vertical', pad=0.01); cbar.set_label(label= 'Rad',fontsize = 20); plt.title('Rad_out', fontsize=20)
        mean_pl, std_pl = np.nanmean(img[:,:,idx_ref], axis=(0,1)), np.nanstd(img[:,:,idx_ref], axis=(0,1)); plt.clim([mean_pl-1*std_pl,mean_pl+2*std_pl]); plt.tight_layout
        
    else:
        fig = plt.figure(figsize=(23,8)); plt.subplot(121); plt.imshow(img[:,:,idx_ref], cmap=cmap); ax = plt.gca(); ax.axis('off'); cbar = plt.colorbar(orientation='vertical', pad=0.01); cbar.set_label(label= 'Rad',fontsize = 20); plt.title('Rad_out', fontsize=20)
        mean_pl, std_pl = np.nanmean(img[:,:,idx_ref], axis=(0,1)), np.nanstd(img[:,:,idx_ref], axis=(0,1)); plt.clim([mean_pl-1*std_pl,mean_pl+2*std_pl]);

        img = mask_threshold (mask_flag_lower, thres_mask_lower, mask_flag_upper, thres_mask_upper, idx_ref, img) #Outputs the same image but with filtered pixels by thresholding
        
        plt.subplot(122); plt.imshow(img[:,:,idx_ref], cmap=cmap); ax = plt.gca(); ax.axis('off'); cbar = plt.colorbar(orientation='vertical', pad=0.01); cbar.set_label(label= 'Rad',fontsize = 20); plt.title('Masked Rad_out', fontsize=20)
        mean_pl, std_pl = np.nanmean(img[:,:,idx_ref], axis=(0,1)), np.nanstd(img[:,:,idx_ref], axis=(0,1)); plt.clim([mean_pl-1*std_pl,mean_pl+2*std_pl]); plt.tight_layout
    
    while plt.fignum_exists(fig.number):
        plt.pause(0.1)  # Esperar un momento antes de verificar nuevamente

    if mission_name == 'EnMAP' or mission_name == 'EMIT':
        return img, wvl_simple, fwhm_simple, sza, vza, idx_ref
    elif mission_name == 'PRISMA' or mission_name == 'GF5':
        return img, wvl_mat, fwhm_mat, sza, vza, idx_ref


### Gas concentration enhancement maps


def gas_processing(img, wvl, fwhm, sza, vza, mission_name, cube, band_names, cont, gas):
    
    h = 0 #km - we ignoring surface height temporary (our default)
    d_h2o = 3 #cm - we ignoring water vapor absorption temporary (our default)
    
    
    if mission_name == 'EnMAP' or mission_name == 'EMIT':
        wvl_simple, fwhm_simple = wvl.copy(), fwhm.copy()
    elif mission_name == 'PRISMA' or mission_name == 'GF5':
        wvl_simple, fwhm_simple = wvl[500].copy(), fwhm[500].copy() #The one in the middle
        
    p_lut_mf = variable_definition.p_lut
        
    window_list = []
    if gas == 'ch4':
        fn_lut_mf = p_lut_mf + 'LUT_ssd_100pm_wvl_2049_2500_nm_ch4.nc' #Extracted with joint combination of ARTS + Libradtran + HITRAN2024 + H2O MT-CKD
        wvl_inf, wvl_sup = 2110, 2450 # As Guanter et al., (2021)
        window = np.where (np.logical_and (wvl_simple >= wvl_inf, wvl_simple <= wvl_sup))[0]; window_list.append(window)
        
    elif gas == 'co2':
        fn_lut_mf = p_lut_mf + 'LUT_ssd_100pm_wvl_1895_2204_nm_co2.nc' #Extracted with joint combination of ARTS + Libradtran + HITRAN2024 + H2O MT-CKD
        wvl_inf, wvl_sup = 1900, 2200 # No need to be accurate. Just to see correlation.
        window = np.where (np.logical_and (wvl_simple >= wvl_inf, wvl_simple <= wvl_sup))[0]; window_list.append(window)
    elif gas == 'h2o':
        fn_lut_mf = p_lut_mf + 'LUT_ssd_100pm_wvl_999_1299_nm_h2o.nc' #Extracted with joint combination of ARTS + Libradtran + HITRAN2024 + H2O MT-CKD
        wvl_inf, wvl_sup = 1000, 1300 #According to Roger et al., (2024), water vapor absorption does not remove the signal
        window = np.where (np.logical_and (wvl_simple >= wvl_inf, wvl_simple <= wvl_sup))[0]; window_list.append(window)
    elif gas == 'c2h4': #Default to SZA = 40º, VZA = 0º, PWV = 3 cm, h = 0 km
        wvl_inf, wvl_sup = 1600, 1700 #nm (where we have data)
        window = np.where (np.logical_and (wvl_simple >= wvl_inf, wvl_simple <= wvl_sup))[0]; window_list.append(window)
    elif gas == 'nh3':
        #Balasus (2026) and Ruzicka (2026) do not agree on the selected windows for nh3. I just use an approximated selection, where I do not really account for water vapor absorption
        fn_lut_mf = p_lut_mf + 'LUT_ssd_100pm_wvl_1399_2500_nm_nh3.nc' #Extracted with joint combination of ARTS + Libradtran + HITRAN2024 + H2O MT-CKD
        wvl_inf, wvl_sup = 2100, 2400 #2300 - This covers the same bands within the main ch4 window. If those are not corrupted by h2o, these neither.
        window = np.where (np.logical_and (wvl_simple >= wvl_inf, wvl_simple <= wvl_sup))[0]; window_list.append(window)
        wvl_inf, wvl_sup = 1900, 2100 #1900
        window = np.where (np.logical_and (wvl_simple >= wvl_inf, wvl_simple <= wvl_sup))[0]; window_list.append(window)
        wvl_inf, wvl_sup = 1400, 1600 #1500
        window = np.where (np.logical_and (wvl_simple >= wvl_inf, wvl_simple <= wvl_sup))[0]; window_list.append(window)
        wvl_inf, wvl_sup = 1400, 2500 #SWIR - Just to see how it performs.
        window = np.where (np.logical_and (wvl_simple >= wvl_inf, wvl_simple <= wvl_sup))[0]; window_list.append(window)
    
    for i_win in range(len(window_list)):
        
        window = window_list[i_win]
        #print(i_win, window[0], window[-1])
        
        if gas == 'c2h4': 
            units = 'ppmm'
            wvl_hr, k_hr = np.load(p_lut_mf + 'wvl_c2h4.npy'), np.load(p_lut_mf + 'k_c2h4.npy')
            _, k_arr = convolution_khr(img, wvl, fwhm, wvl_hr, k_hr, window, mission_name)
        else:
            if gas == 'h2o':
                units = 'cm'
                wvl_hr, rad_hr_arr, delta_x_arr = read_luts_libradtran_h2o(fn_lut_mf, sza, vza)
            else:
                units = 'ppmm'
                wvl_hr, rad_hr_arr, delta_x_arr = read_luts_libradtran(fn_lut_mf, sza, vza, d_h2o, h)
            wvl_ret, k_arr = get_k_libradtran(img, wvl, fwhm, wvl_hr, rad_hr_arr, delta_x_arr, window, mission_name)
            
        ret = mf_retrieval(img, k_arr, window, mission_name, bool_libradtran=True) #ppmm
        if gas == 'nh3' and i_win == 0:
            name_extension = '-2300nm'
        elif gas == 'nh3' and i_win == 1:
            name_extension = '-1900nm'
        elif gas == 'nh3' and i_win == 2:
            name_extension = '-1500nm'
        elif gas == 'nh3' and i_win == 3:
            name_extension = '-SWIR'
        else:
            name_extension = ''
        
        cube = np.dstack((cube, ret)); band_names[f'{gas}-MF{name_extension}({units})'] = f'Band_{cont}'; cont += 1
        cube = np.dstack((cube, medfilt2d(ret, kernel_size = 3))); band_names[f'FILTERED({gas}-MF{name_extension})'] = f'Band_{cont}'; cont += 1
    
    return cube, band_names, cont

def retrieval_maps (img, wvl, fwhm, sza, vza, path_save, mission_name, idx_ref, n):

    #We already have the image ready to be processed
    #Now we need the weighting function (K): it has to be deduced
    #Ideally, I'd like to have a huge LUT for all the gases under different conditions. For that purpose, I'd need to run tons of libradtran scripts. I can do that. So, we would be trusting in LUTs.


    gas_list = ['ch4', 'co2', 'h2o', 'c2h4', 'nh3'] #For now we only do methane
    
    cont = 1
    M,N,B = img.shape; 
    cube = np.zeros((M,N,1)); 
    cube[:,:,0] = img[:,:,idx_ref]; band_names = {'Rad_out': f'Band_{cont}'}; cont += 1
    
    for gas in gas_list:
        print(gas)
        cube, band_names, cont = gas_processing(img, wvl, fwhm, sza, vza, mission_name, cube, band_names, cont, gas)

    md = {'band names': band_names}
    simple_bsq_save_md(cube, path_save, n +'_tool4', md)
    
    return cube



def deltax_rets(p, n, path_save): #Main function to extract the gas concentration enhancement maps

    if n[:4] == 'EMIT':
        mission_name = 'EMIT'
    elif n[:5] == 'ENMAP':
        mission_name = 'EnMAP'
    elif n[:3] == 'PRS':
        mission_name = 'PRISMA'
    elif n[:4] == 'GF5A':
        mission_name = 'GF5'
    else:
        print('Please, enter valid data.')

    print(n)
    
    
    bool_match = True
    match_repeat = False
    
    content = os.listdir(path_save)
    key = n +'_tool4'; len_key = len(key)
    for c in content:
        if c[:len_key] == key:
            bool_match = True
        
    if bool_match:
        clarity_match_bool = False        
        
        while clarity_match_bool == False:
            decision = input('There is already a retrieval from this scene. Do you want to reprocess it? y/n: ')
            if decision == 'n':
                clarity_match_bool = True
            elif decision == 'y':
                match_repeat = True
                clarity_match_bool = True
            else:
                print('Not valid answer. Repeat again.')
          
    if bool_match == False | (bool_match == True & match_repeat == True):
                
        
        bool_mask = False # First we observe Radiance data without masking
        mask_flag_lower, thres_mask_lower, mask_flag_upper, thres_mask_upper = False, 0, False, 0
        img, wvl, fwhm, sza, vza, idx_ref = rad_preview(p, n, bool_mask, mask_flag_lower, thres_mask_lower, mask_flag_upper, thres_mask_upper, mission_name) 
        
        approval = False; bool_mask = True
        while approval == False:
            #Demand masking variables
            mask_flag_lower, thres_mask_lower, mask_flag_upper, thres_mask_upper = masking_asking()
            #print(mask_flag_lower, thres_mask_lower, mask_flag_upper, thres_mask_upper)
            
            if mask_flag_lower == False and mask_flag_upper == False:
                clarity = True; approval = True
            else:
                clarity = False
                img, wvl, fwhm, sza, vza, idx_ref = rad_preview(p, n, bool_mask, mask_flag_lower, thres_mask_lower, mask_flag_upper, thres_mask_upper, mission_name)
                
            while clarity == False:
                decision = input('Do you accept the masking? y/n: ')
                if decision == 'y':
                    approval = True; clarity = True
                elif decision == 'n':
                    clarity = True
                else: 
                    pass   
    
        cube = retrieval_maps (img, wvl, fwhm, sza, vza, path_save, mission_name, idx_ref, n)
        print('Processed!')
        
    else:
        cube = load_hdr_file(path_save, n +'_tool4')
        print('Uploaded!')
    
    return cube, mission_name
