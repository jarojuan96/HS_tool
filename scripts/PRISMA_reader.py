import numpy as np
import h5py
import os

def read_prisma_img(path, filename, swir_flag):
    
    file_name = path + filename + '.he5' 
    
    if swir_flag:
        swir_cube_dat = '/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/SWIR_Cube' ##Va a la metadata y toma los datos del SWIR
        swir_lab = 'Swir'
    else:
        swir_cube_dat = '/HDFEOS/SWATHS/PRS_L1_HCO/Data Fields/VNIR_Cube'
        swir_lab = 'Vnir'

    f = h5py.File(file_name, 'r')

    dset = f[swir_cube_dat]

    ltoa_img = np.flip(np.transpose(dset[:, :, :], axes=[0, 2, 1]), axis=0) 
    dset = f['/KDP_AUX/Cw_' + swir_lab + '_Matrix'] 
    wvl_mat_ini = dset[:, :]
    dset = f['/KDP_AUX/Fwhm_' + swir_lab + '_Matrix'] 
    fwhm_mat_ini = dset[:, :]

    sc_fac = f.attrs['ScaleFactor_' + swir_lab]
    of_fac = f.attrs['Offset_' + swir_lab] 
    vza = 0. #Default since it is not shown in L1 data (it is in L2)
    sza = f.attrs['Sun_zenith_angle'] 

    ltoa_img = ltoa_img / sc_fac - of_fac 

    #Lambda
    wvl_mat_ini = np.flip(wvl_mat_ini, axis=1) 
    li_no0 = np.where(wvl_mat_ini[100, :] > 0)[0] #Select only valid bands
    wvl_mat = np.copy(wvl_mat_ini[:, li_no0])
    wl_center_ini = np.mean(wvl_mat, axis=0)

    #FWHM
    fwhm_mat_ini = np.flip(fwhm_mat_ini, axis=1)
    fwhm_mat = np.copy(fwhm_mat_ini[:, li_no0])
    
    _, _, B_tot = ltoa_img.shape
    #We now make adjustments: obtained after careful inspection
    if swir_flag: 
        if B_tot == len(wl_center_ini):
            ltoa_img = np.flip(ltoa_img, axis=2) 
        else:
            ltoa_img = np.flip(ltoa_img[:, :, :-2], axis=2) 
            
    else:
        if B_tot == len(wl_center_ini):
            ltoa_img = np.flip(ltoa_img, axis=2) 
        else:
            ltoa_img = np.flip(ltoa_img[:, :, 3:], axis=2) #Check it

    _, _, B = ltoa_img.shape
    
    #Just look at this if there is an error of this kind - no apparent corrupt bands
    #idx_arr_discard_fail = np.append(np.arange(0,152), np.arange(153, B)) #Idxs obtained after careful inspection. There was a corrupted band in this position. I think it was corrupt only in the newest versions. Check (how can we identify it without making the retrieval?)
    #ltoa_img = ltoa_img[:, :, idx_arr_discard_fail]       
        
    ltoa_img = np.transpose(ltoa_img, (1,0,2)) #To be aligned with EnMAP and EMIT

    return ltoa_img, wvl_mat, fwhm_mat, vza, sza


def PRISMA_loc_time_gcps(path, filename): 
    

    file_name = path + filename + '.he5' 
    
        
    swir_lat = '/HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Latitude_SWIR'
    swir_lon = '/HDFEOS/SWATHS/PRS_L1_HCO/Geolocation Fields/Longitude_SWIR'
    os.system('h5dump -n 1 ' + file_name + ' > /dev/null 2>&1')

    with h5py.File(file_name, mode='r') as f:

        dset = f[swir_lat]
        lat = np.flip(dset[:,:], axis=0); lat_c = np.mean(lat)
        dset = f[swir_lon]
        lon = np.flip(dset[:,:], axis=0); lon_c = np.mean(lon)


    lines_once = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<Placemarks>', '</Placemarks>']

    M,N = lat.shape[0], lat.shape[1]

    with open(path + 'placemark_' + filename + '.placemark', 'w') as f:

        f.write(lines_once[0])
        f.write('\n')
        f.write(lines_once[1])
        f.write('\n')
        
        cont=0
        for i in range(1,M,98):
            for j in range(1,N,98):
                cont += 1
                lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat[i,j]) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon[i,j]) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(i-0.5) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(j-0.5) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
                f.writelines(lines_loop)
        f.write(lines_once[2])


    time = filename[16:30]

    return time, lat_c, lon_c, lat, lon #We extrac this lat, lon because it already provides the georreference



