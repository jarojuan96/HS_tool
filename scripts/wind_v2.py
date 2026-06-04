#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 17:03:59 2026

@author: javier
"""

## Based on 'wind_v1' from 'https://github.com/jarojuan96/HS_CH4_Emission_Quantification.git'

import numpy as np
import os
from netCDF4 import Dataset
import requests
from scipy import interpolate


def find_nearest_id(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def format_dw_url(year, month, day, hour):
    template = "https://portal.nccs.nasa.gov/datashare/gmao/geos-fp/das/Y{y}/M{m:02d}/D{d:02d}/GEOS.fp.asm.tavg1_2d_slv_Nx.{y}{m:02d}{d:02d}_{h:02d}30.V01.nc4"
    return template.format(y=year, m=month, d=day, h=hour)

def read_string_number(string_number):
    number = 0
    signal = 0
    sign = 1
    if string_number[0] == '-':
        sign = -1
        string_number = string_number[1:]
    for i in range(len(string_number)):
        if signal == 0:
            if string_number[i] != '.':
                number = (10**(len(string_number)-i-1))*int(string_number[i]) + number
            else:
                factor = 10**(len(string_number)-i)
                signal = 1
        elif signal == 1:
            number = (10**(len(string_number)-i))*int(string_number[i]) + number
    number = number*sign/factor
    return number


def get_geos_u10(date:str, hour:int, coords:tuple, pbase): #it stores the file with wind information and then it removes the file
    year, month, day = date.split("-")
    lat_c, lon_c = coords

    dw_url = format_dw_url(year, int(month), int(day), hour)
    
    tmp_path = os.path.join(pbase, "tmp") #Creates a file to temporary store wind data
    os.makedirs(tmp_path, exist_ok=True)
    
    tmp_file = os.path.join(pbase + 'tmp/', "wind.nc")
    r = requests.get(dw_url)
    open(tmp_file, "wb").write(r.content)

    winds = Dataset(tmp_file)

    lat = winds["lat"][:]
    lat_idx = find_nearest_id(lat, lat_c)
    lon = winds["lon"][:]
    lon_idx = find_nearest_id(lon, lon_c)

    u = winds["U10M"][:][0, lat_idx, lon_idx] #horizontal dimension (X-axis)
    v = winds["V10M"][:][0, lat_idx, lon_idx] #vertical dimension (Y-axis)

    winds.close()
    os.remove(tmp_file)

    return u, v

def wind_speed_manual(ts, lat, lon, pbase, psave, name): #YYYYMMDDHHMMSS, lat,lon from metadata

    content_match = os.listdir(psave)
    key = name + '_u_arr.npy'; len_key = len(key)
    bool_match = False
    for c in content_match:
        if c[:len_key] == key:
            bool_match = True
            
    if bool_match == True:
        u_arr_final = np.load(psave + name + '_u_arr.npy')
        u10 = u_arr_final[0]
        print('Wind data uploaded!')

    else:
        ts = str(ts)
        u_array, i = np.zeros((1,9)), 0
        year, month, day, hour, minute = ts[:4], ts[4:6], ts[6:8], int(ts[8:10]), int(ts[10:12])
        
        if hour == 0:
            
            if int(day) < 11:
                day_new = int(day)-1
                day_new = '0' + str(day_new)
            else:
                day_new = int(day)-1
                day_new = str(day_new)
                
            date = year + '-' + month + '-' + day_new 
            hour_new = 23
            
            print('Starting GEOS-FP wind request')
            u_array[i,0], u_array[i,3] = get_geos_u10(date, hour_new, (lat,lon), pbase)
            print('1/3')
            date = year + '-' + month + '-' + day
            u_array[i,1], u_array[i,4] = get_geos_u10(date, hour, (lat,lon), pbase)
            print('2/3')
            u_array[i,2], u_array[i,5] = get_geos_u10(date, hour+1, (lat,lon), pbase)
            print('3/3')
            
        elif hour == 23:
            
            if int(day) < 9:
                day_new = int(day)+1
                day_new = '0' + str(day_new)
            else:
                day_new = int(day)+1
                day_new = str(day_new)
                
            hour_new = 0            
            
            date = year + '-' + month + '-' + day
            print('Starting GEOS-FP wind request')
            u_array[i,0], u_array[i,3] = get_geos_u10(date, hour-1, (lat,lon), pbase)
            print('1/3')
            u_array[i,1], u_array[i,4] = get_geos_u10(date, hour, (lat,lon), pbase)
            print('2/3')
            date = year + '-' + month + '-' + day_new
            u_array[i,2], u_array[i,5] = get_geos_u10(date, hour_new, (lat,lon), pbase)
            print('3/3')
        else:
            
            date = year + '-' + month + '-' + day
            print('Starting GEOS-FP wind request')
            u_array[i,0], u_array[i,3] = get_geos_u10(date, hour-1, (lat,lon), pbase)
            print('1/3')
            u_array[i,1], u_array[i,4] = get_geos_u10(date, hour, (lat,lon), pbase)
            print('2/3')
            u_array[i,2], u_array[i,5] = get_geos_u10(date, hour+1, (lat,lon), pbase)
            print('3/3')
            
        arr_x, arr_y = np.array([hour-0.5, hour+0.5, hour+1.5]), np.array([u_array[i,0], u_array[i,1], u_array[i,2]])
        fu = interpolate.interp1d(arr_x, arr_y, fill_value="extrapolate")
        arr_x, arr_y = np.array([hour-0.5, hour+0.5, hour+1.5]), np.array([u_array[i,3], u_array[i,4], u_array[i,5]])
        fv = interpolate.interp1d(arr_x, arr_y, fill_value="extrapolate")
        final_t = hour + (minute-10)/60
        u_array[i,6], u_array[i,7] = fu(final_t), fv(final_t)
        u_array[i,8] = np.sqrt((u_array[i,6]**2)+(u_array[i,7]**2))
        added_time = (minute-10)/60
        if added_time >= 0:
            u10, ux_hour, uy_hour, ux_add, uy_add = u_array[i,8], u_array[i,1], u_array[i,4], u_array[i,2], u_array[i,5]
        else:
            u10, ux_hour, uy_hour, ux_add, uy_add = u_array[i,8], u_array[i,1], u_array[i,4], u_array[i,0], u_array[i,3]
        u_arr_final = np.array([u10,ux_hour,uy_hour,ux_add,uy_add])
        np.save(psave + name + '_u_arr.npy', u_arr_final)
        
        print('(u_0, v_0) = (' + str(round(u_arr_final[1],2)) + ', ' + str(round(u_arr_final[2],2)) + ') m/s (point 0 for interpolation)')
        print('(u_1, v_1) = (' + str(round(u_arr_final[3],2)) + ', ' + str(round(u_arr_final[4],2)) + ') m/s (point 1 for interpolation)')
        print('u10 = ' + str(round(u_arr_final[0],2)) + ' m/s (used in IME-quantification)')
        
        u10 = u_arr_final[0] #This is the one we are going to use for quantification
    
    return u10




