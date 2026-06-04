import numpy as np
import xml.etree.ElementTree as ET
from scipy.signal import medfilt2d
import spectral.io.envi as envi

def read_ang (tree, lab_ang):
    
    for fact in tree.iter(tag = lab_ang):
        try:
            ang1 = np.float64(fact.find('upper_left').text)
            ang2 = np.float64(fact.find('upper_right').text)
            ang3 = np.float64(fact.find('lower_left').text)
            ang4 = np.float64(fact.find('lower_right').text)
        except:
            pass

    return np.mean([ang1, ang2, ang3, ang4])


def read_xml (file_xml, swir_flg):
    
    tree = ET.parse(file_xml) #read in the XML

    
    if swir_flg:
        lab_str_m = 'swir'    
    else:    
        lab_str_m = 'vnir'    
    
    wl_center = []
    wl_fwhm = []
    gain_arr = []
    offs_arr = []
    for fact in tree.iter(tag = 'bandID'):
        try:
            wvl = np.float64(fact.find('wavelengthCenterOfBand').text)
            fwhm = np.float64(fact.find('FWHMOfBand').text)
            gain = np.float64(fact.find('GainOfBand').text)
            offset = np.float64(fact.find('OffsetOfBand').text)

            wl_center = np.append(wl_center, wvl)
            wl_fwhm = np.append(wl_fwhm, fwhm)
            gain_arr = np.append(gain_arr, gain)
            offs_arr = np.append(offs_arr, offset)
        except:
            pass
    
    
    for fact in tree.iter(tag = lab_str_m + 'ProductQuality'):  
        try:
            num_bd_sp = (np.array(fact.find('numChannelsExpected').text)).astype(int)
        except:
            pass

    num_bd = len(wl_center)
    
    if lab_str_m == 'swir':
        wl_center = wl_center[num_bd-num_bd_sp:]
        wl_fwhm = wl_fwhm[num_bd-num_bd_sp:]
        gain_arr = gain_arr[num_bd-num_bd_sp:]
        offs_arr = offs_arr[num_bd-num_bd_sp:]
    else:   
        wl_center = wl_center[0:num_bd_sp]
        wl_fwhm = wl_fwhm[0:num_bd_sp]
        gain_arr = gain_arr[0:num_bd_sp]
        offs_arr = offs_arr[0:num_bd_sp]

    sza = 90. - read_ang (tree, 'sunElevationAngle')
    saa = read_ang (tree, 'sunAzimuthAngle')
    vaa = read_ang (tree, 'sceneAzimuthAngle')
    vza = read_ang (tree, 'acrossOffNadirAngle')

    for fact in tree.iter(tag = 'specific'):
        try:
            hsf = np.float64(fact.find('meanGroundElevation').text)
        except:
            pass

    return wl_center, wl_fwhm, hsf, sza, saa, vaa, vza, gain_arr, offs_arr 

def import_enmap (file_img_env, file_img_xml, swir_flg, tiff_flg):
    
    #sc_coef = 1.e-2
    sc_coef = 1.e+3
    
    wl_center, wl_fwhm, hsf, sza, saa, vaa, vza, gain_arr, offs_arr = read_xml (file_img_xml, swir_flg)
    
    if tiff_flg:
        import georasters as gr   
        im_tmp = np.transpose(gr.load_tiff(file_img_env), (2, 1, 0))
    else:    
        img = envi.open(file_img_env)
        im_tmp = np.transpose(img.open_memmap(writeable = True), (1, 0, 2))

    ncols, nrows, num_bd = im_tmp.shape

    im_ini = np.zeros([ncols, nrows, num_bd])        
    for bd in range(0, len(wl_center)):        
        im_ini[:, :, bd] = (im_tmp[:, :, bd] * gain_arr[bd] + offs_arr[bd])* sc_coef
        
    im_tmp = 0
    
    return im_ini, wl_center, wl_fwhm, hsf, sza, saa, vaa, vza, ncols, nrows, num_bd




def EnMAP_loc_time_gcps(path_img, name_img):
    
    file_xml = path_img + name_img + '-METADATA.XML'
    
    tree = ET.parse(file_xml)
    a = tree.find('product/navigation/RPC')
    b = tree.find('base/spatialCoverage/boundingPolygon')
    c = tree.find('product/image/swir/dimension')
    lat_off_vnir = np.float64(a[0][2].text)
    lon_off_vnir = np.float64(a[0][3].text)
    lat_off_swir = np.float64(a[187][2].text)
    lon_off_swir = np.float64(a[187][3].text)
    delta_lat = lat_off_vnir - lat_off_swir
    delta_lon = lon_off_vnir - lon_off_swir
    lat_ul = np.float32(b[0][1].text) - delta_lat
    lon_ul = np.float32(b[0][2].text) - delta_lon
    lat_ll = np.float32(b[1][1].text) - delta_lat
    lon_ll = np.float32(b[1][2].text) - delta_lon
    lat_lr = np.float32(b[2][1].text) - delta_lat
    lon_lr = np.float32(b[2][2].text) - delta_lon
    lat_ur = np.float32(b[3][1].text) - delta_lat
    lon_ur = np.float32(b[3][2].text) - delta_lon
    lat_c = np.float32(b[5][1].text) - delta_lat
    lon_c = np.float32(b[5][2].text) - delta_lon
    cols = np.int64(c[0].text); 
    rows = np.int64(c[1].text); 

    with open(path_img + 'placemark_' +  name_img + '.placemark', 'w') as f:

        lines_once = ['<?xml version="1.0" encoding="ISO-8859-1"?>', '<Placemarks>', '</Placemarks>']
        f.write(lines_once[0])
        f.write('\n')
        f.write(lines_once[1])
        f.write('\n')

        cont = 1
        lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_ul) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_ul) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(0) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(0) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
        f.writelines(lines_loop)
        cont = 2
        lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_ll) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_ll) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(0) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(rows) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
        f.writelines(lines_loop)
        cont = 3
        lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_lr) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_lr) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(cols) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(rows) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
        f.writelines(lines_loop)
        cont = 4
        lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_ur) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_ur) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(cols) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(0) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
        f.writelines(lines_loop)
        cont = 5
        lines_loop = ['\t<Placemark name="' + 'gcp_' + str(cont) + '">','\n', '\t\t<LABEL>GCP ' + str(cont) + '</LABEL>', '\n', '\t\t<DESCRIPTION />', '\n', '\t\t<LATITUDE>' + str(lat_c) + '</LATITUDE>', '\n', '\t\t<LONGITUDE>' + str(lon_c) + '</LONGITUDE>', '\n', '\t\t<PIXEL_X>' + str(int(cols/2)) + '</PIXEL_X>', '\n', '\t\t<PIXEL_Y>' + str(int(rows/2)) + '</PIXEL_Y>', '\n', '\t\t<STYLE_CSS>symbol:plus; stroke:#ff8800; stroke-opacity:0.8; stroke-width:1.0</STYLE_CSS>', '\n', '\t</Placemark>', '\n']
        f.writelines(lines_loop)
        f.write(lines_once[2])
        
    time_a, time_b = name_img[29:37], name_img[38:44]
    time = time_a + time_b
    
    lat, lon = None, None

    return time, lat_c, lon_c, lat, lon


def read_enmap_img(path_dat, img_str, swir_flg):
    
    try:
        tiff_flg = True
        file_img_env = path_dat + img_str + '-SPECTRAL_IMAGE_SWIR' + '.TIF'
        file_xml = path_dat + img_str + '-METADATA.XML'
        ltoa_img, wl_center, wl_fwhm, hsf, sza, _, _, vza, ncols, nrows, num_bd = import_enmap (file_img_env, file_xml, swir_flg, tiff_flg)
    except:
        tiff_flg = False #Less common - only shortly after launch
        file_img_env = path_dat + img_str + '-SPECTRAL_IMAGE_SWIR' + '.HDR'
        file_xml = path_dat + img_str + '-METADATA.XML'
        ltoa_img, wl_center, wl_fwhm, hsf, sza, _, _, vza, ncols, nrows, num_bd = import_enmap (file_img_env, file_xml, swir_flg, tiff_flg)

    ltoa_img = np.transpose(ltoa_img, (1,0,2))
    
    return ltoa_img, wl_center, wl_fwhm, hsf, sza, vza


#path_img = '/home1/jroger/Desktop/postdoc/paper_nh3_c2h4/images/images_nh3/Uzbequistan/enmap/ENMAP01-____L1B-DT0000142147_20250706T064239Z_009_V010506_20260319T033150Z/'
#name_img = 'ENMAP01-____L1B-DT0000142147_20250706T064239Z_009_V010506_20260319T033150Z'

#EnMAP_loc_time_gcps(path_img, name_img)

