# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 14:57:15 2018

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
import cmocean

def plot_slices(times,fig_s=(12,5), save_file="autum_stratifications",vmax=None,vmin=None,plot_contour=False, contour_levels=None,bg_color='gray'):
    slices=len(times)
    fig_s=(12,5)
    
    
    fig, ax = plt.subplots(1, slices, sharey=True,figsize=fig_s)
    xlabel=None
    ylabel=None
    show_cb=False
    for fig_slice in range(slices):
        if(fig_slice==slices-1):
            show_cb=True
        start_time=mp.dates.date2num(datetime.datetime.strptime(times[fig_slice][0],"%Y%m%d"))
        end_time=mp.dates.date2num(datetime.datetime.strptime(times[fig_slice][1],"%Y%m%d"))
        plt.axes(ax[fig_slice])
        pfdat.plot_full_data(   new_fig=None,value="temp", save_file="autum_stratifications", \
                                col_map=cmocean.cm.thermal,show_colorbar=show_cb,\
                                tmin=start_time,tmax=end_time,xlabel=xlabel,ylabel=ylabel, \
                                vmax=vmax, vmin=vmin, plot_contour=plot_contour, contour_levels=contour_levels,background_color=bg_color)
        xlabel=""
        ylabel=""
    