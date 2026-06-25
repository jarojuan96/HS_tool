#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:10:15 2026

@author: jroger
"""

import numpy as np
import xml.dom.minidom as md
from osgeo import gdal, osr
import rasterio
from rasterio.transform import xy
import os

from envi_files import read_hdr_bands, load_hdr_file

from EnMAP_reader import EnMAP_loc_time_gcps
from EMIT_reader import EMIT_loc_time_gcps
from PRISMA_reader import PRISMA_loc_time_gcps
from GF_reader import GF_loc_time_gcps



##### Georreferencing and plume source location

def location_and_time(path_folder, name, mission): #Used to create GCPs to georreference data and also to provide timestamp and lat, lon of scene to extract the GEOS-FP wind speed
    
    if mission == 'GF5':
        time, lat_c, lon_c, lat, lon = GF_loc_time_gcps(path_folder, name)
    elif mission == 'EnMAP':
        time, lat_c, lon_c, lat, lon = EnMAP_loc_time_gcps(path_folder, name)
    elif mission == 'PRISMA':
        time, lat_c, lon_c, lat, lon = PRISMA_loc_time_gcps(path_folder, name) 
    elif mission == 'EMIT':
        time, lat_c, lon_c, lat, lon = EMIT_loc_time_gcps(path_folder, name)
        
    return time, lat_c, lon_c, lat, lon

def gcp_from_placemark(placemark_full_path): #Extract GCPs from image, which will later be used to georreference the images

    tree = md.parse(placemark_full_path)
    placemarks = tree.getElementsByTagName("Placemark")
    gcps = []
    for pm in placemarks:
        name = pm.getAttribute("name")
        lat = pm.getElementsByTagName("LATITUDE")[0].childNodes[0].data
        lon = pm.getElementsByTagName("LONGITUDE")[0].childNodes[0].data
        x = pm.getElementsByTagName("PIXEL_X")[0].childNodes[0].data
        y = pm.getElementsByTagName("PIXEL_Y")[0].childNodes[0].data
        gcps.append({"name":name, "lat":lat, "lon":lon, "x":x, "y":y})
        
    return [gdal.GCP(float(p["lon"]), float(p["lat"]), 1, float(p["x"]), float(p["y"])) for p in gcps]


def create_tif_qgis_bands(path_img, name_img, path_ret, name_ret, source_coord): #Ouput bands are georreferenced into Tif files

    array_names = read_hdr_bands(path_ret, name_ret)
    cube = load_hdr_file(path_ret, name_ret)
    M,N,B = cube.shape
    
    if os.path.exists(path_ret + 'georreference_' + name_ret):
        pass
    else:
        os.mkdir(path_ret + 'georreference_' + name_ret)
    
    placemark_full_path = path_img + 'placemark_' + name_img + '.placemark' #stored next to the original L1 data
    list_gcp = gcp_from_placemark(placemark_full_path)
    
    for i in range(B+1):
    
        if i == B: #As a tracer for the pixel source
            fn_tif = path_ret + 'georreference_' + name_ret + '/pixel_source.tif'
            plot = np.ones(cube[:,:,0].shape)
            plot[source_coord[1], source_coord[0]] = 1000000
        else:
            fn_tif = path_ret + 'georreference_' + name_ret + '/' + array_names[i] + '.tif'
            plot = cube[:,:,i].copy()
            
        key_exception = 'GF5'; len_key = len(key_exception)
        if name_ret[:len_key] == key_exception:
            plot = np.flip(np.flip(plot.copy(), axis=1), axis=0)
            N_new = N
            M_new = M
        else:
            N_new = N
            M_new = M
    
        #if i!=0: #Why?
        #    plot[plot<0] = 0    
        data = plot.copy()
    
        driver = gdal.GetDriverByName("GTiff")
        ds = driver.Create(fn_tif , N_new, M_new, 1, gdal.GDT_Float32)    
        band = ds.GetRasterBand(1)
        band.WriteArray(data)
    
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326) #Lat/Lon (deg)
        proj = srs.ExportToWkt()
    
        ds.SetGCPs(list_gcp, proj)
        
        gdal.Warp(fn_tif, ds, dstSRS="EPSG:4326")
        
        ds = None; band = None
    
    return 

def create_tif_qgis_plume(path_img, name_img, path_ret, name_ret, gas, mask): #Plumes are georreferenced into Tif files
    
    array_names = read_hdr_bands(path_ret, name_ret)
    cube = load_hdr_file(path_ret, name_ret)
    M,N,B = cube.shape
    
    if os.path.exists(path_ret + 'georreference_' + name_ret):
        pass
    else:
        os.mkdir(path_ret + 'georreference_' + name_ret)
    
    placemark_full_path = path_img + 'placemark_' + name_img + '.placemark' #stored next to the original L1 data
    list_gcp = gcp_from_placemark(placemark_full_path)
    
    if gas == 'nh3':
        key_match = 'nh3-MF-2300nm(ppmm)'
    elif gas == 'ch4':
        key_match = 'ch4-MF(ppmm)'
    elif gas == 'c2h4':
        key_match = 'c2h4-MF-2300nm(ppmm)'
    
    for i in range(B):
        
        if array_names[i] == key_match:
            
            #print('Match!')
            
            plot = np.ones(cube[:,:,0].shape)*np.nan
            mask = np.float32(mask*1);
            mask[mask==0.] = np.nan
            plot = mask*cube[:,:,i]
    
            fn_tif = path_ret + 'georreference_' + name_ret + f'/{gas}_plume_mask_' + array_names[i] + '.tif'
            
            #if i!=0: #Why?
            #    plot[plot<0] = 0    
            data = plot.copy()
            
            key_exception = 'GF5'; len_key = len(key_exception)
            if name_ret[:len_key] == key_exception:
                data = np.flip(np.flip(data.copy(), axis=1), axis=0)
                N_new = N
                M_new = M
            else:
                N_new = N
                M_new = M
        
            driver = gdal.GetDriverByName("GTiff")
            ds = driver.Create(fn_tif , N_new, M_new, 1, gdal.GDT_Float32)    
            band = ds.GetRasterBand(1)
            band.WriteArray(data)
        
            srs = osr.SpatialReference()
            srs.ImportFromEPSG(4326) #Lat/Lon (deg)
            proj = srs.ExportToWkt()
        
            ds.SetGCPs(list_gcp, proj)
            ds = None; band = None
    
    return 


def coordinates_source(source_coord, path_ret, name_ret): #Use a tracer we made to deduce the plume source location
    
    fn = path_ret + 'georreference_' + name_ret + '/pixel_source.tif'    

    with rasterio.open(fn) as src:
        data = src.read(1)
        
        #print(np.argmax(data))
        idx = np.unravel_index(np.argmax(data), data.shape)
        row, col = idx

    lon_s, lat_s = xy(src.transform, row, col)
    #print(lat_s, lon_s)
    return lat_s, lon_s


def georreference(path_img, name_img, path_ret, name_ret, gas, mask, source_coord): #General function for georreferencing
    
    create_tif_qgis_bands(path_img, name_img, path_ret, name_ret, source_coord)
    create_tif_qgis_plume(path_img, name_img, path_ret, name_ret, gas, mask)
    
    lat, lon = coordinates_source(source_coord, path_ret, name_ret)
    
    return lat, lon

