#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  6 16:10:15 2026

@author: jroger
"""

import numpy as np
from spectral import envi
import re

##### ENVI (.hdr + .bsq) files

def simple_bsq_save_md(img, path, img_name, md): #save ENVI files
    envi.save_image(path + img_name + '.hdr', img, force=True, ext='.bsq',
    interleave='bsq', dtype=np.float32, metadata=md) 
    
def load_hdr_file(path, name): #load ENVI files
    
    img = envi.open(path + name + '.hdr')
    img = img.open_memmap(writeable = True)
    
    return img

def read_hdr_bands(p, n): #load the band names from the ENVI files (retrievals)
    
    with open(p+n+'.hdr', 'r') as f:
        hdr = f.read()

    match = re.search(
        r'band names\s*=\s*\{(.*?)\}',
        hdr,
        re.DOTALL | re.IGNORECASE
    )

    if match:
        band_names = [
            b.strip()
            for b in match.group(1).split(',')
        ]

    return band_names

