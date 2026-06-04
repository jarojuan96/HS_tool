#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 29 17:05:27 2026

@author: jroger
"""

from main import main
import os


def path_name_out (p_gas): #used to process all the images 
    
    p_arr, n_arr, mission_arr, site_arr = [], [], [], []
    
    content_sites = os.listdir(p_gas)
    
    for site in content_sites:
        p_site = p_gas + site + '/'
        content = os.listdir(p_site)
        
        for c in content:
            
            if c == 'emit' or c == 'enmap' or c == 'prisma' or c == 'gf5':
                p = p_site + c + '/'
                
                content_sub = os.listdir(p)
                #print(content_sub)
                
                for c_sub in content_sub:
                    
                    if c == 'enmap':
                        key1 = 'ENMAP01'; len_key1 = len(key1)
                        key2 = 'Z'
                        if c_sub[:len_key1] == key1 and c_sub[-1] == key2:
                            n_arr.append(c_sub)
                            p_arr.append(p + c_sub + '/')
                            mission_arr.append('EnMAP')
                            site_arr.append(site)
                    
                    elif c == 'emit':
                        key1 = 'EMIT_L1B_RAD'; len_key1 = len(key1)
                        if c_sub[:len_key1] == key1:
                            key_skip = '.nc'; len_skip = len(key_skip)
                            n_arr.append(c_sub[:-len_skip])
                            p_arr.append(p)
                            mission_arr.append('EMIT')
                            site_arr.append(site)
                            
                    elif c == 'prisma':
                        key1 = 'PRS_L1_STD'; len_key1 = len(key1)
                        if c_sub[:len_key1] == key1:
                            key_skip = '.he5'; len_skip = len(key_skip)
                            n_arr.append(c_sub[:-len_skip])
                            p_arr.append(p)
                            mission_arr.append('PRISMA')
                            site_arr.append(site)
                            
                    elif c == 'gf5':
                        key1 = 'GF5A_AHSI'; len_key1 = len(key1)
                        if c_sub[:len_key1] == key1:
                            n_arr.append(c_sub)
                            p_arr.append(p + c_sub + '/')
                            mission_arr.append('GF5')
                            site_arr.append(site)
                            
    return p_arr, n_arr, mission_arr, site_arr


p_base = '/home1/jroger/Desktop/postdoc/paper_nh3_c2h4/images/'
gas_arr = ['nh3']

for gas in gas_arr:

    p_gas = f'{p_base}images_{gas}/'
    p_out = f'{p_base}output_{gas}/'
    p_arr, n_arr, mission_arr, site_arr = path_name_out (p_gas)

    for i in range(len(p_arr)):
        
        main(p_arr[i], n_arr[i], p_out, gas, site_arr[i]) 
    


