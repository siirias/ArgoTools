# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 17:13:27 2017

@author: siirias
"""
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np

from plot_depth_n import plot_depth_n
import argohelper as ah
save_file="mbi_depth"
f,(ax1,ax2)=plt.subplots(2,1,sharex=True,figsize=[10,8])
plt.axes(ax1)
line1,=plot_depth_n(new_fig=None,value="salt",target_pressure=180,color='r')
line2,=plot_depth_n(value="salt",target_pressure=150,color='b')
line3,=plot_depth_n(value="salt",target_pressure=120,color='g')
ylims=plt.ylim()
marker_time=mp.dates.date2num(datetime.datetime.strptime("20140821","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)
marker_time=mp.dates.date2num(datetime.datetime.strptime("20150805","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)
plt.ylim(ylims[0],ylims[1])

plt.legend([line1,line2,line3],['180 m','150 m','120 m'])
plt.axes(ax2)

line1,=plot_depth_n(new_fig=None,value="oxygen",target_pressure=180,color='r')
line2,=plot_depth_n(value="oxygen",target_pressure=150,color='b')
line3,=plot_depth_n(value="oxygen",target_pressure=120,color='g')
ylims=plt.ylim()
marker_time=mp.dates.date2num(datetime.datetime.strptime("20140821","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)
marker_time=mp.dates.date2num(datetime.datetime.strptime("20150805","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)
plt.ylim(ylims[0],ylims[1])
#plt.legend([line1,line2],['150 m','120 m'])

if(save_file!=None):
    plt.savefig(save_file+'.png',dpi=300)
    plt.savefig(save_file+'.eps',dpi=300)
#    plt.savefig(save_file+'.jpg',dpi=300)
