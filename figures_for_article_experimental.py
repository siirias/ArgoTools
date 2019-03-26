# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 11:21:35 2016

Plot some figures for the article.
@author: siirias
"""
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np

targe_distance=10
if(len(plt.get_fignums())>0):
    plt.clf()
#get the scripts
runfile('./plot_depth.py', wdir='.')
runfile('./plot_full_data.py', wdir='.')


#plot_full_data(new_fig=True,ref_point=(lats[i],lons[i]), ref_dist=target_distance,time_highlight=date_str ,value="salt", vmin=11,vmax=13,tmin=735600,tmax=735850)
start_time=mp.dates.date2num(datetime.datetime.strptime("20150101","%Y%m%d"))
#end_time=mp.dates.date2num(datetime.datetime.strptime("20150920","%Y%m%d"))
end_time=mp.dates.date2num(datetime.datetime.strptime("20151225","%Y%m%d"))
f_size=(15,10)

print("MBI specific salt and oxygen intersections")
plot_full_data(new_fig=f_size,value="salt", vmin=11,vmax=13,save_file="salt_bottom.png")#,tmin=start_time,tmax=end_time)
plot_full_data(new_fig=f_size,value="salt", vmin=6,vmax=10,save_file="salt_surface.png")#),tmin=start_time,tmax=end_time)

plot_full_data(new_fig=f_size,value="oxygen", vmin=0,vmax=50,save_file="oxygen_bottom.png")#,tmin=start_time,tmax=end_time)
plot_full_data(new_fig=f_size,value="oxygen",save_file="oxygen_full.png")#,tmin=start_time,tmax=end_time)

minbi_start_time=mp.dates.date2num(datetime.datetime.strptime("20160201","%Y%m%d"))
minbi_end_time=mp.dates.date2num(datetime.datetime.strptime("20160801","%Y%m%d"))

plot_full_data(new_fig=f_size,value="oxygen", vmin=0,vmax=50,save_file="minor_inflow_oxygen.png",tmin=minbi_start_time,tmax=minbi_end_time)


plot_full_data(new_fig=f_size,value="temp",save_file="temperature_full.png")#,tmin=start_time,tmax=end_time)
plot_full_data(new_fig=f_size,value="temp", vmin=4,vmax=10,save_file="temperature_bottom.png")#,tmin=start_time,tmax=end_time)

plot_full_data(new_fig=f_size,value="scatter",save_file="scatter_full.png")#,tmin=start_time,tmax=end_time)
plot_full_data(new_fig=f_size,value="scatter",save_file="scatter_small.png",vmin=0.0,vmax=0.001)

"""
print("Full data plots")
plot_full_data(new_fig=True,value="salt")

plot_full_data(new_fig=True,value="temp")

plot_full_data(new_fig=True,value="oxygen")


print("All profiles and elimination")
runfile('./plot_profile_image.py', wdir='.')

print("Routes and area")

runfile('./plot_bathymetry_and_routes.py', wdir='.')


print("MBI Specific")
runfile('./plot_mbi_specific_plots.py', wdir='.')

print("Distance analysis")
plot_full_distance(new_fig=(15,12), ref_point=(57.315,20.04))

runfile('./distance_figure.py', wdir='.')
"""
