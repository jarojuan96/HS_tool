#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:15:11 2026

@author: jroger
"""

import numpy as np
import netCDF4 as nc
from scipy import interpolate
from scipy.interpolate import RegularGridInterpolator


def generate_filter (wvl_M, wvl, wl_resol): #Function for ISRF's convolution
    
    num_wvl_M = wvl_M.shape[0]
    num_wvl = wvl.shape[0]

    s_norm_M = np.zeros((num_wvl_M, num_wvl))
    exp_max = 2.
    exp_min = 2.
    exp_arr = exp_max + (exp_min - exp_max) * np.arange(num_wvl) / (num_wvl-1)
    c_arr = np.power(1./((np.power(2, exp_arr)*np.log(2.))), 1./exp_arr)

    for bd in range(0, num_wvl):

        li1 = np.where (np.logical_and (wvl_M >= (wvl[bd] - 2.* wl_resol[bd]), wvl_M <= (wvl[bd] + 2.* wl_resol[bd])))

        if len(li1[0]) > 0:

            tmp = np.absolute(wvl[bd] - wvl_M[li1])/(wl_resol[bd] *c_arr[bd])
            s = np.exp(-(np.power(tmp, exp_arr[bd])))
            s_norm_M[li1, bd] = s / np.sum(s)
    
    return(s_norm_M)



###-------------------------------
### MODTRAN-based (Luis' LUT - transmittance)


def read_luts_ch4_modtran_trans(file_lut, sza, vza):

    amf = 1. / np.cos(vza * np.pi /180) + 1. / np.cos(sza * np.pi /180) #Assuming geometrical AMF

    t_arr_str, sc_arr_str, mr_arr_str = 't_ch4_arr', 'ch4_sc_arr', 'mr_ch4_arr'

    nc_lut = nc.Dataset(file_lut, 'r', format='NETCDF4') 

    wvl_mod = np.array(nc_lut.variables['wvl_mod']) 

    tmp = nc_lut.variables[t_arr_str] 
    t_arr = np.copy(tmp).T 

    sc_arr = np.array(nc_lut.variables[sc_arr_str])
    tmp = np.array(nc_lut.variables[mr_arr_str])
    mr_arr_all = np.copy(tmp).T

    amf_arr = np.array(nc_lut.variables['amf_arr']) 
    
    f_t = interpolate.interp1d(amf_arr, t_arr, axis = 0)
    f_mr = interpolate.interp1d(amf_arr, mr_arr_all, axis = 0)

    t_arr = f_t(amf)
    mr_arr = f_mr(amf) 

    mr_lim = 2000 #ppb of enhancement (more than enough for what we typically see
    idx_subset = np.where(mr_arr <= mr_lim)[0]

    #print(len(mr_arr[idx_subset]))
    #print(mr_arr[idx_subset])
    
    return wvl_mod, t_arr[idx_subset], sc_arr[idx_subset], mr_arr[idx_subset] #This is set by default for the 'output_Tch4_LUT_AMF_VZA_0_v2.nc' LUT deltaXCH4 sampling



def calc_jac_rad (mr_gas_arr, n_wvl, t_gas_arr, delta_mr_ref):   

    n_pts = len(mr_gas_arr) 
    delta_rad = np.zeros((n_pts, n_wvl)) 
    delta_mr = np.zeros(n_pts) 
    
    for i in range(n_pts):
        
        delta_rad[i, :] = t_gas_arr[i, :] / t_gas_arr[0, :] 
        delta_mr[i] = mr_gas_arr[i] - mr_gas_arr[0]    

    n_pol_jac = 2 
    jac_gas = np.zeros((n_pol_jac+1, n_wvl)) 
    
    for i in range(n_wvl):    
        
        jac_gas[:, i] = np.polyfit(delta_mr, delta_rad[:, i], n_pol_jac)

    mr_poly = np.array([]) 
    for i in range(0, n_pol_jac+1):
        mr_poly = np.append(mr_poly, [delta_mr**i]) 
    
    mr_poly = np.reshape(mr_poly, [n_pol_jac+1, n_pts]) 
    mr_poly = np.flip(mr_poly, axis=0)      


    jac_spec_rad = np.zeros(n_wvl)
    for i in range(n_wvl):    
        jac_spec_rad[i] =  2*jac_gas[0, i] * delta_mr_ref + jac_gas[1, i] 


    return jac_spec_rad


def get_k_ch4_modtran_trans(img, wvl, fwhm, sza, vza, window, fn_lut, mission_name):

    #Default one - Before updating to the new one - Just to make it work (RTM are not specific, we use trans instead of rad, convolution strategy is not accurate)
    wvl_mod, t_gas_arr, _, mr_gas_arr = read_luts_ch4_modtran_trans(fn_lut, sza, vza)
    delta_mr_ref, n_wvl = 1., len(wvl_mod) #ppm, #number of HR-bands
    k_hr = calc_jac_rad (mr_gas_arr, n_wvl, t_gas_arr, delta_mr_ref)
    row, col, B = img.shape

    if mission_name == 'EnMAP' or mission_name == 'EMIT':
        s = generate_filter(wvl_mod, wvl, fwhm)
        k_arr = np.dot(k_hr, s)

    elif mission_name == 'PRISMA' or mission_name == 'GF5':
        k_arr = np.zeros((col, B)) #Non-uniformity in the AT-detector-array
        for i in range(0,col):
            s = generate_filter(wvl_mod, wvl[i], fwhm[i])
            k_arr[i] = np.dot(k_hr, s)
            
    return k_arr



###-------------------------------
#### Libradtran + ARTS

def read_luts_libradtran(file_lut, sza, vza, d_h2o=3, h=0): #We set these defaults due to the extended spectra case

    #We must find the SZA_eff that, together with VZA_eff = 0, is equivalent to the path length (geometrical) from SZA and VZA, assuming: AMF = 1/cos(SZA) + 1/cos(VZA) = 1/cos(SZA_eff) + 1/cos(VZA_eff) = 1/cos(SZA_eff) + 1
    if file_lut == '/home1/jroger/Desktop/postdoc/libRadtran-2.0.6/lut_nc/LUT_nh3_SWIR_v2.nc':
        sza_eff = sza
    else:
        amf = (1 / np.cos(np.deg2rad(sza))) + (1 / np.cos(np.deg2rad(vza)))
        sza_eff = np.rad2deg(np.arccos(1 / (amf - 1)))
        vza = 0

    nc_lut = nc.Dataset(file_lut, 'r', format='NETCDF4') 

    wvl_arr = np.array(nc_lut.variables['wvl_lut'][:]) 
    delta_x_arr = np.array(nc_lut.variables['delta_x_lut'][:])
    h_arr = np.array(nc_lut.variables['h_lut'][:])
    sza_arr = np.array(nc_lut.variables['sza_lut'][:])
    d_h2o_arr = np.array(nc_lut.variables['d_h2o_lut'][:])
    vza_arr = np.array(nc_lut.variables['vza_lut'][:])
    rad_arr = np.array(nc_lut.variables['rad_lut'][:]) # dims: h x sza x d_h2o x vza x delta_x x wvl

    #print(sza_arr, h_arr, d_h2o_arr, delta_x_arr, vza_arr)
    #print('h_arr: ', h_arr, 'km')
    #print('sza_arr: ', sza_arr, 'º')
    #print('d_h2o_arr: ', d_h2o_arr, 'cm')
    #print('vza_arr: ', vza_arr, 'º')
    
    interp = RegularGridInterpolator((h_arr, sza_arr, d_h2o_arr, vza_arr), rad_arr, method='linear') #We set linear as a first approach
    input_values = [h, sza_eff, d_h2o, vza]
    point = [input_values] 
    rad_interp = interp(point); #print(rad_interp.shape)

    rad_interp = rad_interp[0] #delta_x x wvl
    #print(rad_interp.shape)
    return wvl_arr, rad_interp, delta_x_arr

def read_luts_libradtran_h2o(file_lut, sza, vza, h=0): #We set these defaults due to the extended spectra case

   
    amf = (1 / np.cos(np.deg2rad(sza))) + (1 / np.cos(np.deg2rad(vza)))
    sza_eff = np.rad2deg(np.arccos(1 / (amf - 1)))
    vza = 0 #included in sza_eff
    delta_x = 0 #we vary d_h2o instead

    nc_lut = nc.Dataset(file_lut, 'r', format='NETCDF4') 

    wvl_arr = np.array(nc_lut.variables['wvl_lut'][:]) 
    delta_x_arr = np.array(nc_lut.variables['delta_x_lut'][:])
    h_arr = np.array(nc_lut.variables['h_lut'][:])
    sza_arr = np.array(nc_lut.variables['sza_lut'][:])
    d_h2o_arr = np.array(nc_lut.variables['d_h2o_lut'][:])
    vza_arr = np.array(nc_lut.variables['vza_lut'][:])
    rad_arr = np.array(nc_lut.variables['rad_lut'][:]) # dims: h x sza x d_h2o x vza x delta_x x wvl
    rad_arr = np.transpose(rad_arr, (0,1,3,4,2,5))

    #print(sza_arr, h_arr, d_h2o_arr, delta_x_arr, vza_arr)
    #print('h_arr: ', h_arr, 'km')
    #print('sza_arr: ', sza_arr, 'º')
    #print('d_h2o_arr: ', d_h2o_arr, 'cm')
    #print('vza_arr: ', vza_arr, 'º')
    
    interp = RegularGridInterpolator((h_arr, sza_arr, vza_arr, delta_x_arr), rad_arr, method='linear') #We set linear as a first approach
    input_values = [h, sza_eff, vza, delta_x]
    point = [input_values] 
    rad_interp = interp(point); #print(rad_interp.shape)

    rad_interp = rad_interp[0] #d_h2o x wvl
    #print(rad_interp.shape)
    return wvl_arr, rad_interp, d_h2o_arr




def get_k_libradtran_basic(wvl_inst, fwhm_inst, wvl_hr, rad_hr_arr, delta_x_arr, window): #Extract K as in Foote (2021) - more physically consistent

    
    
    wvl, fwhm = wvl_inst[window], fwhm_inst[window] #We get all the bands from the 'valid range', but the boundaries in case there is no enough HR data for convolution at these points
    s = generate_filter(wvl_hr, wvl, fwhm) #Convolution kernel
    
    #print(wvl_inf, wvl_sup)

    rad_conv_arr = []
    for i_delta_x in range(len(delta_x_arr)):

        rad_conv = np.dot(rad_hr_arr[i_delta_x], s)
        rad_conv_arr.append(rad_conv)
        
    rad_conv_arr = np.array(rad_conv_arr)
    
    lograd = np.log(rad_conv_arr.T, out=np.zeros_like(rad_conv_arr.T), where=rad_conv_arr.T > 0); 
    slope, _, _, _ = np.linalg.lstsq(np.stack((np.ones_like(delta_x_arr), delta_x_arr)).T, lograd.T, rcond=None)
    k = slope[1, :]

    return wvl, k

def get_k_libradtran(img, wvl_inst_arr, fwhm_inst_arr, wvl_hr, rad_hr_arr, delta_x_arr, window, mission):
    

    if mission == 'PRISMA' or mission == 'GF5':
        N, B = wvl_inst_arr.shape
        
        k_matrix = []
        wvl_ret = []
        for i_col in range(N):
            wvl, k = get_k_libradtran_basic(wvl_inst_arr[i_col], fwhm_inst_arr[i_col], wvl_hr, rad_hr_arr, delta_x_arr, window)
            k_matrix.append(k); wvl_ret.append(wvl);

    elif mission == 'EnMAP' or mission == 'EMIT' or mission == 'AVIRIS-NG':

        wvl_ret, k_matrix = get_k_libradtran_basic(wvl_inst_arr, fwhm_inst_arr, wvl_hr, rad_hr_arr, delta_x_arr, window)

    return wvl_ret, k_matrix #Different dimensions depending whether it is PRISMA or EnMAP/EMIT


def convolution_khr(img, wvl_inst_arr, fwhm_inst_arr, wvl_hr, k_hr, window, mission):
    
    if mission == 'PRISMA' or mission == 'GF5':
        N, B = wvl_inst_arr.shape
        
        k_matrix = []
        wvl_ret = []
        for i_col in range(N):
            s = generate_filter(wvl_hr, wvl_inst_arr[i_col][window], fwhm_inst_arr[i_col][window])
            k = np.dot(k_hr, s)
            k_matrix.append(k); wvl_ret.append(wvl_inst_arr[i_col][window]);

    elif mission == 'EnMAP' or mission == 'EMIT' or mission == 'AVIRIS-NG':
        
        s = generate_filter(wvl_hr, wvl_inst_arr[window], fwhm_inst_arr[window])
        k_matrix = np.dot(k_hr, s)
        wvl_ret = wvl_inst_arr[window]
        #print(k_matrix)
        
    return wvl_ret, k_matrix #Different dimensions depending whether it is PRISMA or EnMAP/EMIT
    