import netCDF4 as nc
import numpy as np


def read_EMIT(path_img, name_img):
    
    name_obs= name_img[:9] + 'OBS' + name_img[12:]
    file_name = path_img + name_img + '.nc'
    file_obs = path_img + name_obs + '.nc'

    ds = nc.Dataset(file_name)
    ds_obs = nc.Dataset(file_obs)
    
    ltoa = np.array(ds['radiance']) * 10 #W m-2 sr-1 mum-1 (as EnMAP/PRISMA)
    try:
        wvl = np.array(ds.groups['sensor_band_parameters']['wavelengths'])
        fwhm = np.array(ds.groups['sensor_band_parameters']['fwhm'])
    except:
        wvl = np.array(ds.groups['sensor_band_parameters']['radiance_wl'])
        fwhm = np.array(ds.groups['sensor_band_parameters']['radiance_fwhm'])
    vza_arr, sza_arr = np.array(ds_obs['obs'])[:,:,2], np.array(ds_obs['obs'])[:,:,4] #Idxs according to metadata
    vza_m, sza_m = np.mean(vza_arr), np.mean(sza_arr) #Mean VZA and SZA

    ds.close()
    ds_obs.close()
    
    return ltoa, wvl, fwhm, sza_m, vza_m


def EMIT_loc_time_gcps(path_img, name_img): 
    
    file_name = path_img + name_img + '.nc'

    ds = nc.Dataset(file_name)
    
    lon = np.array(ds['location']['lon'])
    lat = np.array(ds['location']['lat'])
    lat_c, lon_c = (np.min(lat) + np.max(lat))/2, (np.min(lon) + np.max(lon))/2 #Mean latitude and longitude

    #name_xml = name_img
    #path_save = path_img + 'xmls/'
    lon = np.array(ds['location']['lon'])
    lat = np.array(ds['location']['lat'])
    with open(path_img + 'placemark_' + name_img + '.placemark', 'w') as f:

        lines_once = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<Placemarks>', '</Placemarks>']
        f.write(lines_once[0])
        f.write('\n')
        f.write(lines_once[1])
        f.write('\n')
        arr_i = np.arange(0,lon.shape[0],100)
        arr_j = np.arange(0,lon.shape[1],100)
        cont = 0
        for i_,i in enumerate(arr_i):
            for j_,j in enumerate(arr_j):
                cont += 1
                lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat[i,j]) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon[i,j]) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(j+1) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(i+1) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
                f.writelines(lines_loop)
        f.write(lines_once[2])
        
    time_a, time_b = name_img[17:25], name_img[26:32]
    time = time_a + time_b
    
    ds.close()
        
    return time, lat_c, lon_c, lat, lon