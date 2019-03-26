# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 11:21:35 2016

Plot some figures for the article.
@author: siirias
"""

import math
import datetime
import matplotlib as mp
import numpy as np
from scipy.io import netcdf


import matplotlib.pyplot as plt
import argohelper as ah
import plot_full_distance as pfdist
import plot_full_data as pfdat
import plot_depth as pdepth
import plot_sal_ox_profiles as psox
import plot_slices as psli
import cmocean
targe_distance=10
#plt.clf()
#get the scripts
#runfile('./plot_depth.py', wdir='.')
#runfile('./plot_full_data.py', wdir='.')


#plot_full_data(new_fig=True,ref_point=(lats[i],lons[i]), ref_dist=target_distance,time_highlight=date_str ,value="salt", vmin=11,vmax=13,tmin=735600,tmax=735850)
start_time=mp.dates.date2num(datetime.datetime.strptime("20150101","%Y%m%d"))
#end_time=mp.dates.date2num(datetime.datetime.strptime("20150920","%Y%m%d"))
end_time=mp.dates.date2num(datetime.datetime.strptime("20151225","%Y%m%d"))

#Plot for the oxygen penetration
contour=True

print("MBI specific salt and oxygen intersections".format( plt.gcf().number))
contour_levels=np.arange(11,13.5,0.1)
pfdat.plot_full_data(new_fig=True,value="salt", col_map=cmocean.cm.haline,vmin=11,vmax=13.5,ymin=100,ymax=220,tmin=start_time,tmax=end_time, save_file="salinity_bottom",plot_contour=contour, contour_levels=contour_levels)

contour_levels=np.arange(0,1.5,0.1)
pfdat.plot_full_data(new_fig=True,value="oxygen",col_map=cmocean.cm.turbid, vmin=0,vmax=1.5,ymin=100,ymax=220,tmin=start_time,tmax=end_time, save_file="oxygen_bottom",plot_contour=contour, contour_levels=contour_levels)
print("MBI specific salt and oxygen intersections (Fig {} and {})".format( plt.gcf().number-1,plt.gcf().number))


contour_levels=np.arange(6,13.5,0.1)
pfdat.plot_full_data(new_fig=True,value="salt", save_file="Salinity_full", col_map=cmocean.cm.haline,plot_contour=contour, contour_levels=contour_levels)

contour_levels=np.arange(0,20.0,0.5)
pfdat.plot_full_data(new_fig=True,value="temp", save_file="Temperature_full", col_map=cmocean.cm.thermal,plot_contour=contour, contour_levels=contour_levels)

contour_levels=np.arange(0,10,0.2)
pfdat.plot_full_data(new_fig=True,value="oxygen",col_map=cmocean.cm.oxy,vmin=0,vmax=10, save_file="Oxygen_full",plot_contour=contour, contour_levels=contour_levels)
print("Full data plots (Figs {},{} and {}".format(plt.gcf().number-2,plt.gcf().number-1,plt.gcf().number))


runfile('./plot_profile_image.py', wdir='.')
print("All profiles and elimination (Fig {})".format(plt.gcf().number))


runfile('./plot_bathymetry_and_routes.py', wdir='.')
print("Routes and area (Fig. {})".format(plt.gcf().number))


runfile('./plot_mbi_specific_plots.py', wdir='.')
print("MBI Specific (Fig {})".format(plt.gcf().number))

psox.plot_sal_ox_profiles('oxygen')
psox.plot_sal_ox_profiles('salt')

pfdist.plot_full_distance(new_fig=(15,6), ref_point=(57.315,20.04))
print("Distance analysis (Fig {})".format(plt.gcf().number))

runfile('./distance_figure.py', wdir='.')

print("Distance figure (Fig {})".format(plt.gcf().number))

runfile('PlotMapWithMeasurements.py')
print("Map and measurements (Fig {})".format(plt.gcf().number))


f_size=(15,6)

#Plot for the oxygen penetration
contour_levels=np.arange(0,2,0.1)
contour=True
oxyg_start_time=mp.dates.date2num(datetime.datetime.strptime("20140821","%Y%m%d"))
oxyg_end_time=mp.dates.date2num(datetime.datetime.strptime("20150101","%Y%m%d"))

pfdat.plot_full_data(   new_fig=f_size,value="oxygen", vmin=0,vmax=2,ymin=50,ymax=200,\
                        tmin=oxyg_start_time,tmax=oxyg_end_time, \
                        save_file="oxygen_penetration",col_map=cmocean.cm.turbid,\
                        plot_contour=contour, contour_levels=contour_levels)
print("Oxygen penetr (Fig {})".format(plt.gcf().number))

f_size=(15,6)
#Some plots for the minor baltic inflow:
minbi_start_time=mp.dates.date2num(datetime.datetime.strptime("20150801","%Y%m%d"))
minbi_end_time=mp.dates.date2num(datetime.datetime.strptime("20160901","%Y%m%d"))
    
    
pfdat.plot_full_data(   new_fig=f_size,value="oxygen",save_file="minor_inflow_oxygen", vmin=0,vmax=1.5,\
                        ymin=100,ymax=220,tmin=minbi_start_time,tmax=minbi_end_time,\
                        col_map=cmocean.cm.turbid,\
                        plot_contour=contour, contour_levels=contour_levels)

times=[["20130901","20131201"],["20141001","20150101"],["20151001","20160101"]]
f_size=(12,5)
contour_levels=np.arange(0,23,0.2)
psli.plot_slices(times,f_size,"autum_stratifications",vmax=None,vmin=None,contour_levels=contour_levels,plot_contour=contour)

runfile('./calibrate_oxygen.py', wdir='.')

ah.give_statistics()

print "FINISHED"