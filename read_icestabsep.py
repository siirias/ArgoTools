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
        points = [line.strip().split(',') for line in file.readlines() if line.strip()]
    points = np.array(points, dtype=float)

    if points.size == 0:
        raise ValueError(f"Polygon file '{file_path}' is empty.")

    if points.shape[1] != 2:
        raise ValueError(f"Polygon file '{file_path}' must contain exactly two columns: latitude, longitude.")

    # The polygon files in this project are stored as (latitude, longitude),
    # but Matplotlib expects (longitude, latitude) for Path.contains_points().
    points = points[:, ::-1]

    # Make sure the polygon is closed by repeating the first point at the end
    if not np.array_equal(points[0], points[-1]):
        points = np.vstack([points, points[0]])

    polygon_path = Path(points)
    return polygon_path

def filter_data_within_polygon(data, polygon_path):
    # Matplotlib Path.contains_points expects (x, y) = (longitude, latitude)
    points = np.column_stack((data['Longitude [degrees_east]'], data['Latitude [degrees_north]']))
    return data[polygon_path.contains_points(points)]

def fill_nan_to_match_length(lists):
    if lists is None or len(lists) == 0:
        return []

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


def resolve_column_name(columns, candidates):
    for name in candidates:
        if name in columns:
            return name
    raise KeyError(f"None of the expected columns were found. Tried: {candidates}")


def parse_ices_dates(values):
    values = values.astype('string').str.strip()
    midnight_next_day = (
        values.str.contains('T24:', regex=False, na=False)
        | values.str.contains('T24Z', regex=False, na=False)
        | values.str.endswith('T24', na=False)
    )
    normalized = values.str.replace('T24:', 'T00:', regex=False)
    normalized = normalized.str.replace('T24Z', 'T00Z', regex=False)
    normalized = normalized.str.replace('T24', 'T00', regex=False)
    dates = pd.to_datetime(normalized, format='mixed', utc=True, errors='coerce')
    return dates + pd.to_timedelta(midnight_next_day.astype('int8'), unit='D')


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
    
    if grouped.empty:
        raise ValueError(
            "No valid ICES records remained after filtering. "
            "Check that the source file contains rows inside the polygon and that the polygon coordinates are valid."
        )

    # Add QCLEVEL and TYPE columns
    grouped['QCLEVEL'] = 'ICES'
    grouped['TYPE'] = 'ICES'

    if grouped['PRES'].empty or grouped['TEMP'].empty or grouped['SAL'].empty:
        raise ValueError(
            "No pressure/temperature/salinity values were found for the filtered data. "
            "This usually means the input file had no usable profile rows after the polygon filter."
        )

    # Fill NaN values in PRES, TEMP, and SAL to match the length of the longest list
    grouped['PRES'] = fill_nan_to_match_length(grouped['PRES'])
    grouped['TEMP'] = fill_nan_to_match_length(grouped['TEMP'])
    grouped['SAL'] = fill_nan_to_match_length(grouped['SAL'])

    return grouped

def read_icestabsep(filename, ofilename, polygon_name, datalines = None):
    if not datalines:
        datalines = [2,-1]
    polygon_path = read_polygon_points(polygon_name)

    columns = pd.read_csv(filename, nrows=0).columns
    longitude_col = resolve_column_name(columns, ['Longitude [degrees_east]'])
    latitude_col = resolve_column_name(columns, ['Latitude [degrees_north]'])
    depth_col = resolve_column_name(columns, ['Depth [m]', 'Depth (ADEPZZ01_ULAA) [m]'])
    pressure_col = resolve_column_name(columns, ['Pressure [dbar]', 'Pressure (PRESPR01_UPDB) [dbar]'])
    temp_col = resolve_column_name(columns, ['Temperature [degC]', 'Temperature (TEMPPR01_UPAA) [degC]'])
    sal_col = resolve_column_name(columns, ['Practical Salinity [dmnless]', 'Salinity (PSALPR01_UUUU) [dmnless]'])

    common_columns = ['Cruise', 'Station', longitude_col, latitude_col,
                      depth_col, pressure_col, temp_col, sal_col]
    if 'yyyy-mm-ddThh:mm:ss.sss' in columns:
        common_columns.append('yyyy-mm-ddThh:mm:ss.sss')
    else:
        common_columns.extend(['Year', 'Month', 'Day', 'Hour', 'Minute'])

    polygon_points = np.loadtxt(polygon_name, delimiter=',', ndmin=2)[:, ::-1]
    polygon_lon = polygon_points[:, 0]
    polygon_lat = polygon_points[:, 1]
    filtered_chunks = []
    n_source_rows = 0
    n_usable_rows = 0
    n_after_filter = 0
    source_lon_min = source_lon_max = source_lat_min = source_lat_max = None

    for chunk in pd.read_csv(filename, usecols=common_columns, chunksize=100000,
                             dtype={'Cruise': 'string', 'Station': 'string'},
                             low_memory=True, on_bad_lines='skip'):
        n_source_rows += len(chunk)
        lon = pd.to_numeric(chunk[longitude_col], errors='coerce')
        lat = pd.to_numeric(chunk[latitude_col], errors='coerce')
        valid_coordinates = lon.notna() & lat.notna()
        if valid_coordinates.any():
            chunk_lon = lon[valid_coordinates]
            chunk_lat = lat[valid_coordinates]
            source_lon_min = chunk_lon.min() if source_lon_min is None else min(source_lon_min, chunk_lon.min())
            source_lon_max = chunk_lon.max() if source_lon_max is None else max(source_lon_max, chunk_lon.max())
            source_lat_min = chunk_lat.min() if source_lat_min is None else min(source_lat_min, chunk_lat.min())
            source_lat_max = chunk_lat.max() if source_lat_max is None else max(source_lat_max, chunk_lat.max())

        depth = pd.to_numeric(chunk[depth_col], errors='coerce')
        pressure = pd.to_numeric(chunk[pressure_col], errors='coerce')
        temp = pd.to_numeric(chunk[temp_col], errors='coerce')
        sal = pd.to_numeric(chunk[sal_col], errors='coerce')
        usable_rows = (depth > 0) | (pressure > 0) | (temp > 0) | (sal > 0)
        chunk = chunk.loc[usable_rows].copy()
        if chunk.empty:
            continue

        n_usable_rows += len(chunk)
        chunk[longitude_col] = lon.loc[chunk.index]
        chunk[latitude_col] = lat.loc[chunk.index]
        inside_polygon = polygon_path.contains_points(
            np.column_stack((chunk[longitude_col].to_numpy(), chunk[latitude_col].to_numpy()))
        )
        chunk = chunk.loc[inside_polygon].copy()
        if chunk.empty:
            continue

        chunk['SOURCEv'] = 'ices_' + chunk['Cruise'].astype(str) + '_' + chunk['Station'].astype(str)
        if 'yyyy-mm-ddThh:mm:ss.sss' in chunk:
            chunk['DATESTR'] = parse_ices_dates(
                chunk['yyyy-mm-ddThh:mm:ss.sss']
            ).dt.strftime('%Y%m%d%H%M%S')
        else:
            chunk['DATESTR'] = (
                chunk['Year'].astype(str) + chunk['Month'].astype(str).str.zfill(2)
                + chunk['Day'].astype(str).str.zfill(2) + chunk['Hour'].astype(str).str.zfill(2)
                + chunk['Minute'].astype(str).str.zfill(2) + '00'
            )
        chunk['DATESv'] = chunk['DATESTR'].astype(float)
        chunk['IDv'] = ('ices_' + chunk['Cruise'].astype(str) + '_' + chunk['Station'].astype(str)
                        + '_' + chunk[longitude_col].astype(str) + '_' + chunk[latitude_col].astype(str))
        filtered_chunks.append(chunk)
        n_after_filter += len(chunk)

    tbl_formatted = pd.concat(filtered_chunks, ignore_index=True) if filtered_chunks else pd.DataFrame()
    n_before_filter = n_usable_rows

    overlap_lon = (source_lon_min is not None and source_lon_min <= polygon_lon.max()
                   and source_lon_max >= polygon_lon.min())
    overlap_lat = (source_lat_min is not None and source_lat_min <= polygon_lat.max()
                   and source_lat_max >= polygon_lat.min())
    print(f"Polygon filter debug: {n_before_filter} usable rows from {n_source_rows} source rows, {n_after_filter} after filter for '{polygon_name}'")
    print(f"Source longitude range: min={source_lon_min}, max={source_lon_max}")
    print(f"Source latitude range: min={source_lat_min}, max={source_lat_max}")
    print(f"Polygon longitude range: min={polygon_lon.min()}, max={polygon_lon.max()}")
    print(f"Polygon latitude range: min={polygon_lat.min()}, max={polygon_lat.max()}")
    print(f"Bounding-box overlap check: lon={overlap_lon}, lat={overlap_lat}, overall={overlap_lon and overlap_lat}")

    if n_after_filter == 0:
        raise ValueError(
            f"Polygon filter removed all rows: {n_before_filter} usable rows -> 0 rows inside polygon '{polygon_name}'. "
            "Check the polygon file and the coordinate ordering used by the dataset."
        )

    tbl_formatted = process_data(tbl_formatted)
    save_to_mat(ofilename, tbl_formatted)    
    
    return tbl_formatted, tbl_formatted

if __name__ == '__main__':
    w_dir = r'/mnt/c/Data/DMQC/UPDATE_test/'  # Work directory
    data_file = 'ICESCTD00-26_v2.csv'
    data_sets = [
        {'source':f'{w_dir}{data_file}','output':f'{w_dir}fmi_ctd_1601.mat', 'polygon':f'{w_dir}polygon_1601.txt'},
        {'source':f'{w_dir}{data_file}','output':f'{w_dir}fmi_ctd_1602.mat', 'polygon':f'{w_dir}polygon_1602.txt'},
        {'source':f'{w_dir}{data_file}','output':f'{w_dir}fmi_ctd_1501.mat', 'polygon':f'{w_dir}polygon_1501.txt'},
        {'source':f'{w_dir}{data_file}','output':f'{w_dir}fmi_ctd_1502.mat', 'polygon':f'{w_dir}polygon_1502.txt'}
        ]
    for d in data_sets:
        tbl, tbl_formatted = read_icestabsep(d['source'],d['output'],d['polygon'])
        print(f"{d['source']} -> {d['output']}: original rows={tbl.shape[0]}, filtered rows={tbl_formatted.shape[0]}")
