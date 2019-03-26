# -*- coding: utf-8 -*-
"""
Created on Thu May 07 16:45:33 2015

@author: siirias
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pylab as p
import re
import os
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from netCDF4 import Dataset

def haversine_dist(lon, lat, lon0, lat0):
    '''Haversine distance of lon and lat from a single point lon0 and lat0.
    Application from http://stackoverflow.com/questions/19413259/efficient-way-to-calculate-distance-matrix-given-latitude-and-longitude-data-in
    '''

    earth_radius = 6371.0 * 1000
    lons = np.deg2rad(lon0)
    lats = np.deg2rad(lat0)
    lon = np.deg2rad(lon)
    lat = np.deg2rad(lat)

    #dlat, dlon = (lats - lat, lons - lon)
    lat_dif = (lats / 2 - lat / 2)
    lon_dif = (lons / 2 - lon / 2)

    #a = np.sin(dlat/2.0)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2.0)**2
    lat_dif=np.sin(lat_dif)
    lon_dif=np.sin(lon_dif)

    lat_dif=np.power(lat_dif, 2)
    lon_dif=np.power(lon_dif, 2)

    lon_dif *= (np.cos(lats) * np.cos(lat))
    dif = lon_dif + lat_dif

    dif=np.arctan2(np.power(dif, .5), np.power(1 - dif, .5))
    dif = dif * (2 * earth_radius)

    return dif  # in metres
        
arrow_cmap = mpl.colors.LinearSegmentedColormap('my_colormap', \
            {'red':((0.0,1.0,1.0),(0.3,0.1,0.1),(1.0,0.0,0.0)), 
           'green':((0.0,0.0,0.0),(0.3,0.5,0.5),(1.0,0.0,0.0)),
           'blue':((0.0,0.0,0.0),(0.3,0.0,0.0),(1.0,0.0,0.0)) } \
           ,256)
#execfile("../../pack_argo_data.py")
#dataA1=np.loadtxt('profilesAPE1.dat')
#dataA2=np.loadtxt('profilesAPE2.dat')
plot_type='point_number'  #point_number, data
press_breaker=90.0 #the depth taht cuts higher, and lower layer.
p_switch='higher' # 'higher', 'lower','none'
#p_switch='lower' # 'higher', 'lower','none'
#,DATE (YYYY/MM/DD HH:MI:SS),LATITUDE (degree_north),LONGITUDE (degree_east),PLATFORM,PRES LEVEL1(decibar),PSAL LEVEL1(psu),QC,TEMP LEVEL1(degree_Celsius),time difference,short td,latp,lonp,Distance (km),Distance_sn (km),Direction,Speed (m/s),shallow,Speed_sn (m/s),cycle mean pressure,floating,no profile

#NewData=np.loadtxt('Apex_traj_dataset_computed_variables.csv',delimiter=',',skiprows=1,usecols=(0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
NewData=np.genfromtxt('Apex_traj_dataset_computed_variables.csv',delimiter=',',skip_header=1)
data_sets=[NewData]
c_lon=3
c_lat=2
c_plat=4
c_dt=9
c_shortdt_flag=10
c_dist=13
c_speed=16
c_floating=20
c_noprof=21
c_press=19
lon2m=52202
lat2m=111194
grid_step_lon=0.05*4
grid_step_lat=0.025*4
if(p_switch=='lower'):
    grid_step_lon=0.05*2
    grid_step_lat=0.025*2

#Let's figure aut the delta movements:
recorded_speeds=[]
move_vecs=[]; #Here comes teh move-vectors between every profile
for data in data_sets:
#    d_t=np.array(data[:,6])
#    for i in range(d_t.size-1):
#        d_t[i]=d_t[i+1]-d_t[i]
    
    for i in range(data[1:,0].size-1):
        d_lon=data[i,c_lon]-data[i-1,c_lon]
        d_lat=data[i,c_lat]-data[i-1,c_lat]
        if(abs(d_lon)+abs(d_lat)>0.001): #this is an actual move
                if(data[i,c_shortdt_flag]!=1 and data[i,c_floating]!=1 and data[i,c_noprof]!=1): #if time is long enough to qualify
                    if(i>0 and data[i,c_plat]!=data[i-1,c_plat]):
                        print "vaihtui!",speed
                    else:
                        if((p_switch=='higher' and data[i,c_press]<=press_breaker) or \
                            (p_switch=='lower' and data[i,c_press]>=press_breaker) or \
                            (p_switch=='none')):
                            d_lon=d_lon/(data[i,c_dt])
                            d_lat=d_lat/(data[i,c_dt])
                            #d_lon is lon_deg/s d_lat is lat_deg/s
        #                    recorded_speeds.append(100*60*np.sqrt((d_lon*lon2m)**2+(d_lat*lat2m)**2))
        #                    speed=data[i,c_dist]/data[i,c_dt]
                            speed=data[i,c_speed]
    #                        if(speed>0.002):
    #                            print speed
                            recorded_speeds.append(speed)
                            tmp=[data[i,c_lon],data[i,c_lat],d_lon,d_lat]
                            move_vecs.append(tmp)
#                else:
#                    print data[i,c_shortdt_flag],data[i,c_floating],data[i,c_noprof]

#fig = plt.figure()
#plt.plot(data[:,6]-data[0,6],d_depth[:],'.-')
#plt.plot(data[:,6]-data[0,6],right_depth[:],'r.-')
#plt.plot(data[:,6]-data[0,6],'b.-')
#plt.plot(data[:,6],data[:,0],'b.-')
#plt.plot(data[:,6],-1*data[:,0],'g.-')
#plt.plot(data[:,6],-1*data[:,8],'r.-')
#plt.plot(d_t[:],'r.-')
#plt.show()

brdr_lon=1.5;
brdr_lat=0.25;
latmin=data[:,c_lat].min();
latmax=data[:,c_lat].max();
lonmin=data[:,c_lon].min();
lonmax=data[:,c_lon].max();
for data in data_sets:
    latmin=np.array([data[:,c_lat].min(), latmin]).min();
    latmax=np.array([data[:,c_lat].max(), latmax]).max();
    lonmin=np.array([data[:,c_lon].min(), lonmin]).min();
    lonmax=np.array([data[:,c_lon].max(), lonmax]).max();
latmax+=brdr_lat
latmin-=brdr_lat
lonmax+=brdr_lon
lonmin-=brdr_lon
#latmax=np.array([data[:,5].max(),data[:,5].max()]).max()+brdr_lat;
#lonmin=np.array([data[:,4].min(),data[:,4].min()]).min()-brdr_lon;
#lonmax=np.array([data[:,4].max(),data[:,4].max()]).max()+brdr_lon;

#Let's build a gridded version of the move vectors.
gridded_vecs=[]
for vec in move_vecs:
    new_lon=np.round((vec[0]-lonmin)/grid_step_lon)*grid_step_lon+lonmin
    new_lat=np.round((vec[1]-latmin)/grid_step_lat)*grid_step_lat+latmin
    old_vec=None;
    for gvec in gridded_vecs:
        #Check here if this gridpoint allready has a vector
        if(abs(gvec[0]-new_lon)+abs(gvec[1]-new_lat)<0.01):
            gvec[2]+=vec[2]
            gvec[3]+=vec[3]
            gvec[4]+=1 #This counts how many vectors contribute to the final
            old_vec=gvec;
    if old_vec==None:
        gridded_vecs.append([new_lon,new_lat,vec[2],vec[3],1])


plt.clf()
m = Basemap(projection='merc',\
     resolution='h', llcrnrlat=latmin, llcrnrlon=lonmin, urcrnrlat=latmax, urcrnrlon=lonmax)
#m = Basemap(projection='merc',\
#     resolution='l', llcrnrlat=-80, llcrnrlon=-170, urcrnrlat=80, urcrnrlon=170)
#Basemap.drawcoastlines(m)
#m.drawparallels(np.arange(58,61,1),labels=[1,0,0,0])
#m.drawmeridians(np.arange(22,31,2),labels=[0,0,0,1])
#m.drawparallels(np.arange(latmin,latmax,0.5),labels=[1,0,0,0])
#m.drawmeridians(np.arange(lonmin,lonmax,1.0),labels=[0,0,0,1])
#m.fillcontinents(color='coral',lake_color='aqua')

url = 'http://ferret.pmel.noaa.gov/thredds/dodsC/data/PMEL/etopo5.nc'
etopodata = Dataset('iowtopo2_rev03.nc')

topoin = etopodata.variables['Z_TOPO'][:]
lons = etopodata.variables['XT_I'][:]
lats = etopodata.variables['YT_J'][:]
# shift data so lons go from -180 to 180 instead of 20 to 380.
#topoin,lons = shiftgrid(180.,topoin,lons,start=False)
m.drawparallels(np.arange(np.floor(latmin),np.ceil(latmax),0.5),labels=[1,0,0,0],linewidth=0.0)
m.drawmeridians(np.arange(np.floor(lonmin),np.ceil(lonmax),1),labels=[0,0,0,1],linewidth=0.0)
m.fillcontinents(color='#cccccc',lake_color='#ccccff')
Basemap.drawcoastlines(m)

nx = int((m.xmax-m.xmin)/5000.)+1; 
ny = int((m.ymax-m.ymin)/5000.)+1
topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
im = m.imshow(topodat,cm.GMT_haxby)
if (plot_type=='data'):
    plt.colorbar(im,boundaries=range(-150,0,10))
#m.bluemarble()
#p.clf();
x,y = m(data[:,4],data[:,5])
#plt.plot(x,y,'ob--')
speeds=[]
"""
for vec in gridded_vecs:
    asx,asy=m(vec[0],vec[1])
    length=2/vec[4] #vec[4] is the amount of vectors summed in each gridpoint
#    aex,aey=m(vec[0]+vec[2]*length*lat2m,vec[1]+vec[3]*length*lon2m)
#    aex,aey=m(vec[0]+0.02*lon2m,vec[1]+0.02*lat2m)
#    plt.arrow(asx,asy,aex-asx,aey-asy, head_width=3000, head_length=5000, fc='r', ec='r')
    plt.arrow(asx,asy,1e6*vec[2]*lon2m/vec[4],1e6*vec[3]*lat2m/vec[4], head_width=3000, head_length=5000, fc='r', ec='r')
    speeds.append(haversine_dist(asx,asy,asx+vec[2]/vec[4],asy+vec[3]/vec[4]))

#plot the legend-hack
px=lonmin+0.8*(lonmax-lonmin)
py=latmin+0.8*(latmax-latmin)
x,y=m(px,py)
dx,dy=m(px+0.001,py)
print dx,dy
print x,y
plt.arrow(x,y,dx-x,dy-y, head_width=3000, head_length=5000, fc='b', ec='k')
plt.text(x,y-30000,r"$10 \frac{cm}{s}$",fontsize=20)    
"""
tstx=[]
tsty=[]
tstu=[]
tstv=[]
tstc=[]
for vec in gridded_vecs:
    if(vec[4]<2):
        coltxt=1
    else:
        coltxt=0
    if(vec[4]>1):  #eliminate vectors with too little values
        print vec[4]
        tstc.append(np.array([vec[4],20]).min());
        tstx.append(vec[0])
        tsty.append(vec[1])
#        tstu.append(vec[2]*lon2m*100*60/(vec[4])) #cm/min
#        tstv.append(vec[3]*lat2m*100*60/(vec[4])) #cm/min
        tstu.append(vec[2]*lon2m*100/(vec[4])) #cm/s
        tstv.append(vec[3]*lat2m*100/(vec[4])) #cm/s
               
"""
for vec in move_vecs:
    tstx.append(vec[0])
    tsty.append(vec[1])
    tstu.append(vec[2]*lon2m*100*60) #cm/min
    tstv.append(vec[3]*lat2m*100*60) #cm/min
    tstc.append(1);
"""

#Draw the grid for gridded vecs
aval=0.2
min_x=np.min(tstx)-grid_step_lon*0.5
max_x=np.max(tstx)+grid_step_lon*0.5
min_y=np.min(tsty)-grid_step_lat*0.5
max_y=np.max(tsty)+grid_step_lat*0.5
for x_grid in np.arange(min_x,max_x+grid_step_lon*0.5,grid_step_lon):
    tmpx,tmpy=m([x_grid, x_grid],[min_y,max_y])
    plt.plot(tmpx,tmpy,'k',alpha=aval)
for y_grid in np.arange(min_y,max_y+grid_step_lat*0.5,grid_step_lat):
    tmpx,tmpy=m([min_x,max_x],[y_grid,y_grid])
    plt.plot(tmpx,tmpy,'k',alpha=aval)

tstx,tsty=m(tstx,tsty)
#cmap=plt.cm.hot;
cmap=arrow_cmap;
Q=plt.quiver(tstx,tsty,tstu,tstv,tstc,width=0.0020,headwidth=5,scale=60, color="k",cmap=cmap)
plt.quiverkey(Q,0.9,0.8,5,r'$4 \frac{cm}{s}$',fontproperties={'size':20})
if (plot_type=='point_number'):
    mn=1
    mx=21
    md=1
    cbar=plt.colorbar(Q,cmap=cmap,boundaries=range(mn,mx,md))
    cbar.set_ticks(range(mx))
#    cbar.set_ticklabels([mn,md,mx])
    name_base='gridded_argo_surface_amounts'
    if(p_switch=='lower'):
        name_base='gridded_argo_bottom_amounts'
plt.savefig(name_base+'.png' ,facecolor='w',dpi=300)
plt.savefig(name_base+'.ps' ,facecolor='w',dpi=300)
plt.savefig(name_base+'.eps' ,facecolor='w',dpi=300)
print 'average speed for the vectors:', np.average(recorded_speeds)

plt.show()
