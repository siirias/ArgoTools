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

line1,=plot_depth_n(new_fig=(15,10),value="salt",target_pressure=150,color='b')
line2,=plot_depth_n(value="salt",target_pressure=120,color='g')
marker_time=mp.dates.date2num(datetime.datetime.strptime("20140821","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)
marker_time=mp.dates.date2num(datetime.datetime.strptime("20150805","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)

plt.legend([line1,line2],['150 m','120 m'])
#plt.title('Salinity')

fig_name='plotted_depth_salt_new'
plt.savefig('%s.png' % (fig_name),dpi=300)
plt.savefig('%s.eps' % (fig_name),dpi=300)

line1,=plot_depth_n(new_fig=(15,10),value="oxygen",target_pressure=150,color='b')
line2,=plot_depth_n(value="oxygen",target_pressure=120,color='g')
marker_time=mp.dates.date2num(datetime.datetime.strptime("20140821","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)
marker_time=mp.dates.date2num(datetime.datetime.strptime("20150805","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)

plt.legend([line1,line2],['150 m','120 m'])
#plt.title('Oxygen')
fig_name='plotted_depth_oxygen_new'
plt.savefig('%s.png' % (fig_name),dpi=300)
plt.savefig('%s.eps' % (fig_name),dpi=300)

