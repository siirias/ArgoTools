# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 17:41:47 2023
@author: siirias
"""
import pandas as pd
import numpy as np
from io import StringIO
from matplotlib.path import Path
import scipy.io


def save_to_mat(file_path, data):
    # Convert DataFrame to a dictionary
    mat_data = {col: data[col].values for col in data.columns}

    # Save the dictionary to a .mat file
    scipy.io.savemat(file_path, mat_data, do_compression=True)
    
    
    
def read_polygon_points(file_path):
    with open(file_path, 'r') as file:
        points = [line.strip().split(',') for line in file.readlines()]
    points = np.array(points, dtype=float)
    # Make sure the polygon is closed by repeating the first point at the end
    if not np.array_equal(points[0], points[-1]):
        points = np.vstack([points, points[0]])

    polygon_path = Path(points)        
    return polygon_path

def filter_data_within_polygon(data, polygon_path):
    # Assuming 'LONG' and 'LAT' are columns in your DataFrame
    points = np.column_stack((data['Latitude [degrees_north]'], data['Longitude [degrees_east]'] ))
    return data[polygon_path.contains_points(points)]

def fill_nan_to_match_length(lists):
    max_length = max(len(lst) for lst in lists)
    return [lst + [float('nan')] * (max_length - len(lst)) for lst in lists]

# Define a function to strip quotation marks
def reformat_csv(filename):
    # Read the file into a list of strings
    with open(filename, 'r') as file:
        lines = file.readlines()
    # Strip quotation marks from each line
    lines = [line.strip('"\n') for line in lines]
    # Convert the list of strings to a file-like object
    data = StringIO('\n'.join(lines))
    return data

def process_data(tbl):
    # Group by 'IDv' and aggregate data
    try: #old dataformat
        grouped = tbl.groupby('IDv').agg({
            'Longitude [degrees_east]': 'first',
            'Latitude [degrees_north]': 'first',
            'DATESv': 'first',
            'Pressure [dbar]': list,
            'Temperature [degC]': list,
            'Practical Salinity [dmnless]': list,
            'SOURCEv': 'first'
        }).reset_index()
    
        # Rename columns as needed
        grouped.rename(columns={
            'Longitude [degrees_east]': 'LONG',
            'Latitude [degrees_north]': 'LAT',
            'DATESv': 'DATES',
            'Pressure [dbar]': 'PRES',
            'Temperature [degC]': 'TEMP',
            'Practical Salinity [dmnless]': 'SAL',
            'SOURCEv': 'SOURCE'
        }, inplace=True)
    except KeyError: #then try new format
        grouped = tbl.groupby('IDv').agg({
            'Longitude [degrees_east]': 'first',
            'Latitude [degrees_north]': 'first',
            'DATESv': 'first',
            'Pressure (PRESPR01_UPDB) [dbar]': list,
            'Temperature (TEMPPR01_UPAA) [degC]': list,
            'Salinity (PSALPR01_UUUU) [dmnless]': list,
            'SOURCEv': 'first'
        }).reset_index()
    
        # Rename columns as needed
        grouped.rename(columns={
            'Longitude [degrees_east]': 'LONG',
            'Latitude [degrees_north]': 'LAT',
            'DATESv': 'DATES',
            'Pressure (PRESPR01_UPDB) [dbar]': 'PRES',
            'Temperature (TEMPPR01_UPAA) [degC]': 'TEMP',
            'Salinity (PSALPR01_UUUU) [dmnless]': 'SAL',
            'SOURCEv': 'SOURCE'
        }, inplace=True)
    
    # Add QCLEVEL and TYPE columns
    grouped['QCLEVEL'] = 'ICES'
    grouped['TYPE'] = 'ICES'

    # Fill NaN values in PRES, TEMP, and SAL to match the length of the longest list
    grouped['PRES'] = fill_nan_to_match_length(grouped['PRES'])
    grouped['TEMP'] = fill_nan_to_match_length(grouped['TEMP'])
    grouped['SAL'] = fill_nan_to_match_length(grouped['SAL'])

    return grouped

def read_icestabsep(filename, ofilename, polygon_name, datalines = None):
    if not datalines:
        datalines = [2,-1]
    
    tbl = pd.read_csv(reformat_csv(filename), sep=None, engine='python', comment='/') #sep=None should work with tabs or commas
    #tbl = pd.read_csv(reformat_csv(filename))
#    tbl = pd.read_csv(reformat_csv(filename), delimiter = '\\t', engine = 'python')
    try: #old format               
        excl = tbl[(tbl['QV:ODV:Depth [m]'] > 0) | 
                   (tbl['QV:ODV:Pressure [dbar]'] > 0) | 
                   (tbl['QV:ODV:Temperature [degC]'] > 0)| 
                   (tbl['QV:ODV:Practical Salinity [dmnless]'] > 0) 
                    ].index
    except: #new format
        excl = tbl[(tbl['QV:ODV:Depth (ADEPZZ01_ULAA) [m]'] > 0) | 
                   (tbl['QV:ODV:Pressure (PRESPR01_UPDB) [dbar]'] > 0) | 
                   (tbl['QV:ODV:Temperature (TEMPPR01_UPAA) [degC]'] > 0)| 
                   (tbl['QV:ODV:Salinity (PSALPR01_UUUU) [dmnless]'] > 0) 
                    ].index
    
    tbl.drop(excl, inplace=True)

    # Create new columns
    tbl['SOURCEv'] = 'ices_' + tbl['Cruise'].astype(str) + '_' + tbl['Station'].astype(str)
    try:
        # Old ICES format
        tbl['DATESTR'] = (
            tbl['Year'].astype(str)
            + tbl['Month'].astype(str).str.zfill(2)
            + tbl['Day'].astype(str).str.zfill(2)
            + tbl['Hour'].astype(str).str.zfill(2)
            + tbl['Minute'].astype(str).str.zfill(2)
            + '00'
        )
    
    except KeyError:
        # New ICES format
        tbl['DATESTR'] = (
            pd.to_datetime(tbl['yyyy-mm-ddThh:mm:ss.sss'],format='mixed',utc=True)
          .dt.strftime('%Y%m%d%H%M%S'))
                
    
    tbl['DATESv'] = tbl['DATESTR'].astype(float)

    tbl['IDv'] = 'ices_' + tbl['Cruise'].astype(str) + '_' + tbl['Station'].astype(str) + '_' + \
                 tbl['Longitude [degrees_east]'].astype(str) + '_' + \
                 tbl['Latitude [degrees_north]'].astype(str)
    
    polygon_path = read_polygon_points(polygon_name)
    tbl_formatted = filter_data_within_polygon(tbl, polygon_path)
    tbl_formatted = process_data(tbl_formatted)
    save_to_mat(ofilename, tbl_formatted)    
    
    return  tbl, tbl_formatted

w_dir = r'C:\Data\DMQC\UPDATE_test\\'  # Work directory
data_file = 'ICESCTD00-26_bottle.csv'
data_sets = [
    {'source':f'{w_dir}{data_file}','output':f'{w_dir}fmi_ctd_1601.mat', 'polygon':f'{w_dir}polygon_1601.txt'},
    {'source':f'{w_dir}{data_file}','output':f'{w_dir}fmi_ctd_1602.mat', 'polygon':f'{w_dir}polygon_1602.txt'},
    {'source':f'{w_dir}{data_file}','output':f'{w_dir}fmi_ctd_1501.mat', 'polygon':f'{w_dir}polygon_1501.txt'},
    {'source':f'{w_dir}{data_file}','output':f'{w_dir}fmi_ctd_1502.mat', 'polygon':f'{w_dir}polygon_1502.txt'}
    ]
for d in data_sets:
    tbl, tbl_formatted = read_icestabsep(d['source'],d['output'],d['polygon'])
    print(tbl)
    print(tbl.shape, tbl_formatted.shape)