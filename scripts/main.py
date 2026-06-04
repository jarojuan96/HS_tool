#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:09:17 2026

@author: jroger
"""

from hs_tool4_deltaxgas import deltax_rets #1) Get gas concentration enhancement maps (ppmm)
from quant_func_v2 import emission_quantification #2) Quantification: plume delineation + wind speed extraction + IME quantification
from store_results import excel_info #3) Store results into a .csv file


#import warnings
#warnings.filterwarnings("ignore")

def main (p, n, psave, gas, site):
    
    print(site)
    
    cube, mission = deltax_rets(p, n, psave) #The retrieval is done as long as it is not already done or we don't want to reprocess
    
    if gas == 'ch4':
        dxgas = cube[:,:,1] #ch4-MF(ppmm)
        rad_ref = cube[:,:,0] #Rad-ref
        Q, err_Q, u10, err_u10, lat_s, lon_s, ts, bool_det = emission_quantification (dxgas, rad_ref, mission, gas, p, n, psave) #Plume detection - integrated in select_plume function within the 'quant_func_v2'
                                                                              #If no plume is found. Then, we stop the process.
                                                                              #If it is found, we extract wind direction if not already stored                                                                    
        if bool_det:
            fields = ['Site', 'Mission', 'Timestamp\nYYYYMMDDhhmmss', 'Source-lat(º)', 'Source-lon(º)', 'u10 (m/s)', 'err(u10)', 'Q (kg/h)', 'err(Q)']
            info = [site, mission, ts, round(lat_s,4), round(lon_s,4), round(u10,2), round(err_u10,2), round(Q,2), round(err_Q,2)]
    
    elif gas == 'nh3':
        dxgas = cube[:,:,9] #nh3-MF-2300nm(ppmm)
        rad_ref = cube[:,:,0] #Rad-ref
        Q_1, err_Q_1, Q_2, err_Q_2, u10, err_u10, lat_s, lon_s, ts, bool_det = emission_quantification (dxgas, rad_ref, mission, gas, p, n, psave) #Q_1 assumes tau = inf and Q_2 tau = 1h
        
        if bool_det:
            fields = ['Site', 'Mission', 'Timestamp\nYYYYMMDDhhmmss', 'Source-lat(º)', 'Source-lon(º)', 'u10 (m/s)', 'err(u10)', 'Q_tau=inf (kg/h)', 'err(Q_tau=inf)', 'Q_tau=1h (kg/h)', 'err(Q_tau=1h)']
            info = [site, mission, ts, round(lat_s,4), round(lon_s,4), round(u10,2), round(err_u10,2), round(Q_1,2), round(err_Q_1,2), round(Q_2,2), round(err_Q_2,2)]
    
    if bool_det:
        excel_info(fields, info, psave, gas)
    
    return

site = 'Iraq' #name of the location from the acquisition
gas = 'nh3' #gas of interest (potential plume detection). Potential gases so far: ch4 and nh3
p = '/home1/jroger/Desktop/postdoc/paper_nh3_c2h4/images/images_nh3/Iraq/emit/' #path to the L1 file
n = 'EMIT_L1B_RAD_001_20250821T075748_2523305_016' #name of the L1 file without the extension
psave = '/home1/jroger/Desktop/postdoc/paper_nh3_c2h4/images/output_nh3/' #where all the outputs are saved

main (p, n, psave, gas, site)

