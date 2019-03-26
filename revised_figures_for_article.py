# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 10:57:36 2018

@author: siirias
"""

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

#print the routes
runfile('./plot_bathymetry_and_routes.py', wdir='.')

contour=True
bg_col='white'
#Print the main data-figures
save_file="main_dataset"
f,(ax1,ax2,ax3)=plt.subplots(3,1,sharex=True,figsize=[10,10])

contour_levels=np.arange(6,13.5,0.1)
plt.axes(ax1)
pfdat.plot_full_data(new_fig=None,value="salt", col_map=cmocean.cm.haline,plot_contour=contour, contour_levels=contour_levels,background_color=bg_col)

contour_levels=np.arange(0,20.0,0.5)
plt.axes(ax2)
pfdat.plot_full_data(new_fig=None,value="temp",  col_map=cmocean.cm.thermal,plot_contour=contour, contour_levels=contour_levels,background_color=bg_col)

contour_levels=np.arange(0,10,0.2)
plt.axes(ax3)
pfdat.plot_full_data(new_fig=None,value="oxygen",col_map=cmocean.cm.oxy,vmin=0,vmax=10, plot_contour=contour, contour_levels=contour_levels,background_color=bg_col)

if(save_file!=None):
    plt.savefig(save_file+'.png',dpi=300)
    plt.savefig(save_file+'.eps',dpi=300)
#    plt.savefig(save_file+'.jpg',dpi=300)
#"""

save_file="distance_from_deployment"
pfdist.plot_full_distance(new_fig=(15,6), ref_point=(57.315,20.04))
print("Distance analysis (Fig {})".format(plt.gcf().number))
if(save_file!=None):
    plt.savefig(save_file+'.png',dpi=300)
    plt.savefig(save_file+'.eps',dpi=300)
#    plt.savefig(save_file+'.jpg',dpi=300)
    
    
#plot_full_data(new_fig=True,ref_point=(lats[i],lons[i]), ref_dist=target_distance,time_highlight=date_str ,value="salt", vmin=11,vmax=13,tmin=735600,tmax=735850)
start_time=mp.dates.date2num(datetime.datetime.strptime("20150101","%Y%m%d"))
#end_time=mp.dates.date2num(datetime.datetime.strptime("20150920","%Y%m%d"))
end_time=mp.dates.date2num(datetime.datetime.strptime("20160901","%Y%m%d"))

#Plot for the oxygen penetration
contour=True
f,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=[10,8])

print("MBI specific salt and oxygen intersections".format( plt.gcf().number))
contour_levels=np.arange(11,13.5,0.1)
plt.axes(ax1)
pfdat.plot_full_data(new_fig=None,value="salt", col_map=cmocean.cm.haline,vmin=11,vmax=13.5,ymin=100,ymax=220,tmin=start_time,tmax=end_time, plot_contour=contour, contour_levels=contour_levels,background_color=bg_col)

contour_levels=np.arange(0,1.5,0.1)
plt.axes(ax2)
pfdat.plot_full_data(new_fig=None,value="oxygen",col_map=cmocean.cm.turbid, vmin=0,vmax=1.5,ymin=100,ymax=220,tmin=start_time,tmax=end_time, plot_contour=contour, contour_levels=contour_levels,background_color=bg_col)
print("MBI specific salt and oxygen intersections (Fig {} and {})".format( plt.gcf().number-1,plt.gcf().number))
save_file="mbi_specific_sal_ox"
if(save_file!=None):
    plt.savefig(save_file+'.png',dpi=300)
    plt.savefig(save_file+'.eps',dpi=300)
#    plt.savefig(save_file+'.jpg',dpi=300)





runfile('./plot_mbi_specific_plots_revised.py', wdir='.')

times=[["20130901","20131201"],["20141001","20150101"],["20151001","20160101"]]
f_size=(12,5)
contour_levels=np.arange(0,23,0.2)
psli.plot_slices(times,f_size,"autum_stratifications",vmax=None,vmin=None,contour_levels=contour_levels,plot_contour=contour,bg_color=bg_col)

runfile('PlotMapWithMeasurements.py')
save_file="CTD_Argo_places"
if(save_file!=None):
    plt.savefig(save_file+'.png',dpi=300)
    plt.savefig(save_file+'.eps',dpi=300)
 #   plt.savefig(save_file+'.jpg',dpi=300)

runfile('statistics_from_argos.py')
save_file="profileamounts_Gotland"
if(save_file!=None):
    plt.savefig(save_file+'.png',dpi=300)
    plt.savefig(save_file+'.eps',dpi=300)


#"""
 
print "Finished!, \n But remember to run the calibrate_oxygen.py, to get the omparison image."
print "nvm. let's just run it"
runfile("calibrate_oxygen.py")

print "Now, we are Finished!"