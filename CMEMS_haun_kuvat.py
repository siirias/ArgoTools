# -*- coding: utf-8 -*-
"""
Created on Wed Oct 25 11:52:12 2017

@author: siirias
"""

plot_full_data(value="oxygen",col_map=cmocean.cm.oxy,plot_contour=True,vmin=0,vmax=10, contour_levels=np.arange(0,10,0.5))


runfile('D:/ArgoData/plot_bathymetry_and_routes_simple.py', wdir='D:/ArgoData')


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