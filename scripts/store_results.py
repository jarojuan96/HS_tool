#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:09:17 2026

@author: jroger
"""

import csv
import os

def excel_info(fields, info, psave, gas): #Creates a csv file with all the information. Emissions for a given gas are recorded in a specific csv file
    
    name_excel = f'log-{gas}-emissions'
    fn = psave + name_excel + '.csv'
    
    
    file_exists = os.path.exists(fn)
    with open(fn, "a", newline="") as f:
        writer = csv.writer(f)
    
        # write header only if file is new or empty
        if not file_exists or os.path.getsize(fn) == 0:
            writer.writerow(fields)

        writer.writerow(info)
        
    return