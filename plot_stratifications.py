# -*- coding: utf-8 -*-
"""
Created on Wed Nov 22 15:28:48 2017

@author: siirias
"""
import datetime
runfile('D:/ArgoData/plot_depth_n.py', wdir='D:/ArgoData')
tmin=datetime.datetime(2014,4,1)
tmax=datetime.datetime(2014,8,1)
tdiff=datetime.timedelta(years=1)
depths=[10,25,50]
c=['b','r','g']
#full image
plt.figure(figsize=(15,10))
for i in range(len(depths)):
    print c[i],depths[i]
    l,=plot_depth_n(value="temp",target_pressure=depths[i],color=c[i]+'*-')
    #l,=plot_depth_n(value="temp",target_pressure=depths[i],color=c[i]+'*-',timeshift=tdiff)


plt.figure(figsize=(15,10))
for i in range(len(depths)):
    print c[i],depths[i]
    l,=plot_depth_n(value="temp",target_pressure=depths[i],color=c[i]+'*-',tmin=tmin,tmax=tmax)
#line2,=plot_depth_n(value="temp",target_pressure=depths[1],color='g')


marker_time=mp.dates.date2num(datetime.datetime.strptime("20140821","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)
marker_time=mp.dates.date2num(datetime.datetime.strptime("20150805","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)

plt.legend([line1,line2],[str(depths[0])+' m',str(depths[1])+' m'])