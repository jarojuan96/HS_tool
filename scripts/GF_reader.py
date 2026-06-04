import numpy as np
import os
import xml.etree.ElementTree as ET
import georasters as gr 


def import_ahsi (file_img_gtif, file_cal, file_srf, file_xml, wl_1, wl_2):

      
    try:
        data = np.genfromtxt(file_srf, skip_header=0, delimiter=',')
        wl_center_ini = np.array(data[:, 0])
        wl_fwhm_ini = np.array(data[:, 1])
    except:
        data = np.genfromtxt(file_srf, skip_header=0, delimiter=' ')
        wl_center_ini = np.array(data[:, 0])
        wl_fwhm_ini = np.array(data[:, 1])

    try:
        data = np.genfromtxt(file_cal, skip_header=0, delimiter=',')
        #data = np.genfromtxt(file_cal, skip_header=0)
        if np.isnan(data[0]):
            data = np.loadtxt(file_cal)
            cal_coef = np.copy(data[:, 0])
        else:
            cal_coef_tmp = np.array(data)
            cal_coef = np.copy(cal_coef_tmp[:, 0])
    except:        
        #cal_coef = np.genfromtxt(file_cal, skip_header=0)
        data = np.genfromtxt(file_cal, skip_header=0, delimiter=',')
        cal_coef_tmp = np.array(data)
        try:
            cal_coef = np.copy(cal_coef_tmp[:, 0])
        except:
            cal_coef = np.copy(cal_coef_tmp)
   
    filename, file_extension = os.path.splitext(file_xml)
    head_tail = os.path.split(filename) 
    str_tmp = head_tail[1]
    if str_tmp[0:4] == 'ZY1E':
        cal_coef_new = np.resize(cal_coef, len(wl_center_ini))
        cal_coef = np.copy(cal_coef_new)
    
    img_tmp = np.asarray(gr.load_tiff(file_img_gtif))
    #if str_tmp[0:4] == 'GF5_' or str_tmp[0:4] == 'GF5A':
    #    img_open = np.transpose(img_tmp, (1, 2, 0))
    #else:        
    #   img_open = np.transpose(img_tmp, (1, 2, 0))
    img_open = np.transpose(img_tmp, (1, 2, 0))
    img_tmp = 0
    
    [nrows, ncols, num_bd] = img_open.shape
    
    li1 = np.where (np.logical_and (wl_center_ini >= wl_1, wl_center_ini <= wl_2))
    
    wl_center = wl_center_ini[li1]
    wl_fwhm = wl_fwhm_ini[li1]
    cal_coef_col = np.squeeze(cal_coef[li1])
    toa_img = np.squeeze(img_open[:, :, li1])
    num_bd = len(wl_center)
    
    img_open = 0
    
    ltoa_img = np.zeros((ncols, nrows, num_bd))
    for ind_c in range(0, ncols):           
        for ind_r in range(0, nrows):
            ltoa_img[ind_c, ind_r, :] = toa_img[ind_r, ind_c, :] * cal_coef_col 


    tree = ET.parse(file_xml) #read in the XML
    for fact in tree.iter(tag = 'ProductMetaData'):
        vza = np.float16(fact.find('SatelliteZenith').text)
        sza = np.float16(fact.find('SolarZenith').text)

    #print('GF5 (VZA, SZA): ', vza, sza)
    """
    if vza >= 10. or vza <= 0.:
        vza = 0.
    if sza <= 10. or sza >= 60.:
        sza = 45.
    """
    wvl_mat = np.tile(wl_center, (ncols, 1))
    fwhm_mat = np.tile(wl_fwhm, (ncols, 1))
    toa_img = 0
    ltoa_img = np.transpose(ltoa_img, (1,0,2))
    #print('nrows x ncols x nbands = ' + str(nrows) + 'x' + str(ncols) + 'x' + str(num_bd))
    
    return np.float32(ltoa_img), np.float32(wvl_mat), np.float32(fwhm_mat), np.float32(sza), np.float32(vza)


def GF_loc_time_gcps(path_folder, name): #it provides information about the acquisition timestamp and location
    
    path_dat = path_folder[:-(len(name)+1)];
    
    try:
        file_xml = path_dat + name + '/' + name + '.xml'
        filename, file_extension = os.path.splitext(file_xml)
        tree = ET.parse(file_xml)
    except:

        file_xml = path_dat + name + '/' + name + '.xml'
        filename, file_extension = os.path.splitext(file_xml)
        tree = ET.parse(file_xml)
        
    for fact in tree.iter(tag = 'ProductMetaData'):
        time = fact.find('StartTime').text
        rows = np.int64(fact.find('LinesInPixels').text)
        cols = np.int64(fact.find('SamplesInPixels').text)
        lat_c = np.float32(fact.find('CenterLatitude').text)
        lon_c = np.float32(fact.find('CenterLongitude').text)
        lat_tl = np.float32(fact.find('TopLeftLatitude').text)
        lon_tl = np.float32(fact.find('TopLeftLongitude').text)
        lat_tr = np.float32(fact.find('TopRightLatitude').text)
        lon_tr = np.float32(fact.find('TopRightLongitude').text)
        lat_br = np.float32(fact.find('BottomRightLatitude').text)
        lon_br = np.float32(fact.find('BottomRightLongitude').text)
        lat_bl = np.float32(fact.find('BottomLeftLatitude').text)
        lon_bl = np.float32(fact.find('BottomLeftLongitude').text)
        
        
        with open(path_folder + 'placemark_' + name + '.placemark', 'w') as f:

            lines_once = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<Placemarks>', '</Placemarks>']
            f.write(lines_once[0])
            f.write('\n')
            f.write(lines_once[1])
            f.write('\n')
            
            cont = 1
            lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_tl) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_tl) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(0) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(0) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
            f.writelines(lines_loop)
            cont = 2
            lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_bl) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_bl) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(0) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(rows) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
            f.writelines(lines_loop)
            cont = 3
            lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_br) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_br) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(cols) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(rows) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
            f.writelines(lines_loop)
            cont = 4
            lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_tr) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_tr) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(cols) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(0) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
            f.writelines(lines_loop)
            cont = 5
            lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_c) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_c) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(int(cols/2)) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(int(rows/2)) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
            f.writelines(lines_loop)
            f.write(lines_once[2])
        
    YYYY,MM,DD,hh,mm,ss = time[0:4], time[5:7], time[8:10], time[11:13], time[14:16], time[17:19]
    time = YYYY+MM+DD+hh+mm+ss
    lat, lon = None, None
            
    return time, lat_c, lon_c, lat, lon


def read_gf_img(path_folder, name):

    path_dat = path_folder[:-(len(name)+1)]; #print('PATH_DAT', path_dat)
    mission = name[:4]
    
    try:
        file_img_env = path_dat + name + '/' + name + '_SW.tif'
        file_cal = path_dat + name + '/' + mission + '_AHSI_RadCal_SWIR.raw'
        file_srf = path_dat + name + '/' + mission + '_AHSI_Spectralresponse_SWIR.raw'
        file_xml = path_dat + name + '/' + name + '.xml'
        img, wvl_mat, fwhm_arr, sza, vza = import_ahsi (file_img_env, file_cal, file_srf, file_xml, 1008, 2450)
    except:
        
        try:
           file_img_env = path_dat + name + '_SW.tif'
           file_cal = path_dat + mission + '_AHSI_RadCal_SWIR.raw'
           file_srf = path_dat + '/' + mission + '_AHSI_Spectralresponse_SWIR.raw'
           file_xml = path_dat + '/' + name + '.xml'
           img, wvl_mat, fwhm_arr, sza, vza = import_ahsi (file_img_env, file_cal, file_srf, file_xml, 1008, 2450) 
        except:
            file_srf = path_dat + name + '/' + mission +'_AHSI_SWIR_Spectralresponse.raw'
            file_img_env = path_dat + name + '/' + name + '_SW.tiff'
            file_cal = path_dat + name + '/' + mission + '_AHSI_SWIR_RadCal.raw'
            file_xml = path_dat + name + '/' + name + '.xml'
            img, wvl_mat, fwhm_arr, sza, vza = import_ahsi (file_img_env, file_cal, file_srf, file_xml, 1008, 2450)
            

    return img, wvl_mat, fwhm_arr, sza, vza
