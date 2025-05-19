# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 13:27:13 2023

@author: siirias

This script reads a netCDF file containing bathymetry data for the Baltic Sea and creates an interactive map
using the Plotly library in Python. The map displays the bathymetry data as markers, with each marker indicating
the depth at a specific point in the sea. The user can interact with the map to zoom in/out and hover over markers
to view the depth information. This script requires the netCDF4, numpy, and plotly libraries to be installed.
"""

import netCDF4 as nc
import numpy as np
import plotly.graph_objs as go
import plotly.io as pio

data_dir = 'C:/data/ArgoData/'
out_dir = 'c:/data/argodata/figures/'

topo_file = 'iowtopo2_rev03.nc'
topo_variable = 'Z_NEAR'
lat_variable = 'YT_J'
lon_variable = 'XT_I'
landmask_variable = 'LANDMASK'
topo_missing_value = -1.0e34

# Load the netCDF file
dataset = nc.Dataset(data_dir + topo_file)

# Extract the bathymetry data
depth = dataset.variables[topo_variable][:]
# Define the latitude and longitude coordinates
latitudes = dataset.variables[lat_variable][:]
longitudes = dataset.variables[lon_variable][:]
# Load the landmask data
landmask = dataset.variables[landmask_variable][:]

depth[landmask>0.5] = np.nan


# Define the center point and zoom level of the map
center_lat, center_lon = latitudes.mean(), longitudes.mean()
zoom_level = 7

# Define the colorscale for the heatmap
colorscale = [
    [0, 'blue'],
    [0.5, 'cyan'],
    [0.6, 'lime'],
    [0.7, 'yellow'],
    [0.8, 'orange'],
    [1, 'red']
]

# Create the heatmap trace
heatmap_trace = go.Heatmap(
    z=depth,
    x=longitudes,
    y=latitudes,
    zauto=False,
#    zmin=np.nanmin(depth),
#    zmax=np.nanmax(depth),
    zmin = -440,
    zmax = 0,
    colorscale=colorscale,
    name='Depth',
    showscale=True,
    opacity=0.6,
)


# Create the layout for the map
layout = go.Layout(
    title='Baltic Sea Bathymetry',
    autosize=True,
    hovermode='closest',
    xaxis=dict(
        title='Longitude',
    ),
    yaxis=dict(
        title='Latitude',
    ),
    mapbox=dict(
        center=dict(
            lat=center_lat,
            lon=center_lon,
        ),
        zoom=zoom_level,
        style='stamen-terrain',
    ),
)

# Create the figure object
fig = go.Figure([heatmap_trace], layout)

# # Add the coastlines
# fig.add_trace(
#     go.Choroplethmapbox(
#         geojson='https://raw.githubusercontent.com/deldersveld/topojson/master/countries/finland/finland-counties.json',
#         locations=[],
#         z=[],
#         text=[],
#         featureidkey='properties.NAME_2',
#         marker_line_width=0,
#         showscale=False,
#     )
# )
fig.data = fig.data[::-1] # hack to reverse plotting order of things.
# Save the plot as an interactive HTML file
pio.write_html(fig, file=out_dir + 'baltic_sea_bathymetry.html', auto_open=False)