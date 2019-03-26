# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 17:13:27 2017

@author: siirias
"""

line1,=plot_depth(new_fig=(15,10),value="salt",target_pressure=150,color='b')
line2,=plot_depth(value="salt",target_pressure=120,color='g')
marker_time=mp.dates.date2num(datetime.datetime.strptime("20140821","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)
marker_time=mp.dates.date2num(datetime.datetime.strptime("20150805","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)


plt.legend([line1,line2],['150 m','120 m'])
