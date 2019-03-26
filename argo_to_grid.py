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
def argo_to_grid(p_switch='lower',p_gludge_switch='none'):

    #p_switch='none' # 'higher', 'lower','none'
    #p_gludge_switch='none' # picks the given cell-size, unless none. 'higher', 'lower','none'

    """        
    arrow_cmap = mpl.colors.LinearSegmentedColormap('my_colormap', \
                {'red':((0.0,1.0,1.0),(0.3,0.1,0.1),(1.0,0.0,0.0)), 
               'green':((0.0,0.0,0.0),(0.3,0.5,0.5),(1.0,0.0,0.0)),
               'blue':((0.0,0.0,0.0),(0.3,0.0,0.0),(1.0,0.0,0.0)) } \
               ,256)
    """
    arrow_cmap = mpl.colors.LinearSegmentedColormap('my_colormap', \
                {'red':((0.0,0.5,0.5),(0.3,0.5,0.5),(1.0,1.0,1.0)), 
               'green':((0.0,0.0,0.0),(0.3,0.5,0.5),(1.0,1.0,1.0)),
               'blue':((0.0,0.0,0.0),(0.3,0.2,0.2),(1.0,0.8,0.8)) } \
               ,256)
    
    
    #execfile("../../pack_argo_data.py")
    #dataA1=np.loadtxt('profilesAPE1.dat')
    #dataA2=np.loadtxt('profilesAPE2.dat')
    plot_type='point_number'  #point_number, data
    press_breaker=90.0 #the depth taht cuts higher, and lower layer.
    #,DATE (YYYY/MM/DD HH:MI:SS),LATITUDE (degree_north),LONGITUDE (degree_east),PLATFORM,PRES LEVEL1(decibar),PSAL LEVEL1(psu),QC,TEMP LEVEL1(degree_Celsius),time difference,short td,latp,lonp,Distance (km),Distance_sn (km),Direction,Speed (m/s),shallow,Speed_sn (m/s),cycle mean pressure,floating,no profile
    
    #NewData=np.loadtxt('Apex_traj_dataset_computed_variables.csv',delimiter=',',skiprows=1,usecols=(0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19))
    #NewData=np.genfromtxt('Apex_traj_dataset_computed_variables.csv',delimiter=',',skip_header=1)
    NewData=pd.read_csv('complete_Argo_dataset_GoB_All.csv')
    data_sets=[NewData]
    """
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
    """
    cn_lon='lon'
    cn_lat='lat'
    cn_plat='WMO'
    cn_dt='Cyclet'
    cn_dist='Distance (km)'
    cn_speed='Velocity'
    cn_press='Ave drift press'
    cn_date='Date'
    cn_time='Time'
    min_members_in_group=3  #How many vectors must be in a cell to show it
    lon2m=52202
    lat2m=111194
    grid_step_lon=0.05*4
    grid_step_lat=0.025*4
    if((p_switch=='lower' and p_gludge_switch!='higher' )or p_gludge_switch=='lower'):
        grid_step_lon=0.05*2
        grid_step_lat=0.025*2
    
    #Let's figure aut the delta movements:
    recorded_speeds=[]
    move_vecs=[]; #Here comes the move-vectors between every profile
    speed=0
    for data in data_sets:
        
        for i in range(1,len(data)):
            d_lon=data[cn_lon][i]-data[cn_lon][i-1]
            d_lat=data[cn_lat][i]-data[cn_lat][i-1]
            datestr=data[cn_date][i]+' '+data[cn_time][i]
            if(abs(d_lon)+abs(d_lat)>0.001): #this is an actual move
                if(i>1 and data[cn_plat][i]!=data[cn_plat][i-1]):
                    print "vaihtui!",speed
                else:
                    if((p_switch=='higher' and data[cn_press][i]<=press_breaker) or \
                        (p_switch=='lower' and data[cn_press][i]>=press_breaker) or \
                        (p_switch=='none')):
                        d_lon=d_lon/(data[cn_dt][i])
                        d_lat=d_lat/(data[cn_dt][i])
                        speed=data[cn_speed][i]
                        recorded_speeds.append(speed)
                        tmp=[data[cn_lon][i],data[cn_lat][i],d_lon,d_lat,datestr]
                        move_vecs.append(tmp)
    
    brdr_lon=1.5;
    brdr_lat=0.25; #0.25
    latmin=data[cn_lat].min();
    latmax=data[cn_lat].max();
    lonmin=data[cn_lon].min();
    lonmax=data[cn_lon].max();
    for data in data_sets:
        latmin=np.array([data[cn_lat].min(), latmin]).min();
        latmax=np.array([data[cn_lat].max(), latmax]).max();
        lonmin=np.array([data[cn_lon].min(), lonmin]).min();
        lonmax=np.array([data[cn_lon].max(), lonmax]).max();
    latmax+=brdr_lat
    latmin-=brdr_lat
    lonmax+=brdr_lon
    lonmin-=brdr_lon
    
    
    #Let's build a gridded version of the move vectors.
    gridded_vecs=[]
    for vec in move_vecs:
        new_lon=np.round((vec[0]-lonmin)/grid_step_lon)*grid_step_lon+lonmin
        new_lat=np.round((vec[1]-latmin)/grid_step_lat)*grid_step_lat+latmin
        old_vec=None;
        for gvec in gridded_vecs:
            #Check here if this gridpoint already has a vector
            if(abs(gvec[0]-new_lon)+abs(gvec[1]-new_lat)<0.01):
                gvec[2]+=vec[2]
                gvec[3]+=vec[3]
                gvec[4]+=1 #This counts how many vectors contribute to the final
                old_vec=gvec;
        if old_vec==None:
            gridded_vecs.append([new_lon,new_lat,vec[2],vec[3],1,len(gridded_vecs)+1])


    #let's build a list of original vectors, stacked based on their cell:
    orig_vec_data=[]
    first=True
    for i in range(len(move_vecs)):
        vec=move_vecs[i]
        new_lon=np.round((vec[0]-lonmin)/grid_step_lon)*grid_step_lon+lonmin
        new_lat=np.round((vec[1]-latmin)/grid_step_lat)*grid_step_lat+latmin
        for j in range(len(gridded_vecs)):
            gvec=gridded_vecs[j]
            if(abs(gvec[0]-new_lon)+abs(gvec[1]-new_lat)<0.01):
                ucomp=(vec[2]*lon2m*100) #cm/s
                vcomp=(vec[3]*lat2m*100) #cm/s
                newline=[gvec[5],vec[1],vec[0],vcomp,ucomp,vec[4]]
                orig_vec_data.append(newline)

    plt.clf()
    m = Basemap(projection='merc',\
         resolution='h', llcrnrlat=latmin, llcrnrlon=lonmin, urcrnrlat=latmax, urcrnrlon=lonmax)
    
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

    #Draw the grid for gridded vecs
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
        if(vec[4]>=min_members_in_group):  #eliminate vectors with too little values
            tstc.append(np.array([vec[4],20]).min());
            tstx.append(vec[0])
            tsty.append(vec[1])
            tstu.append(vec[2]*lon2m*100/(vec[4])) #cm/s
            tstv.append(vec[3]*lat2m*100/(vec[4])) #cm/s


    aval=1.0
    grid_color='#8080a0'
    linewidth='0.25'
    min_x=np.min(tstx)-grid_step_lon*0.5
    max_x=np.max(tstx)+grid_step_lon*0.5
    min_y=np.min(tsty)-grid_step_lat*0.5
    max_y=np.max(tsty)+grid_step_lat*0.5
    for x_grid in np.arange(min_x,max_x+grid_step_lon*0.5,grid_step_lon):
        tmpx,tmpy=m([x_grid, x_grid],[min_y,max_y])
        plt.plot(tmpx,tmpy,grid_color,alpha=aval,linewidth=linewidth)
    for y_grid in np.arange(min_y,max_y+grid_step_lat*0.5,grid_step_lat):
        tmpx,tmpy=m([min_x,max_x],[y_grid,y_grid])
        plt.plot(tmpx,tmpy,grid_color,alpha=aval,linewidth=linewidth)

    
    nx = int((m.xmax-m.xmin)/5000.)+1; 
    ny = int((m.ymax-m.ymin)/5000.)+1
    topodat = m.transform_scalar(topoin,lons,lats,nx,ny)
    im = m.imshow(topodat,cm.GMT_haxby)
    if (plot_type=='data'):
        plt.colorbar(im,boundaries=range(-150,0,10))
    #m.bluemarble()
    #p.clf();
    #x,y = m(data[cn_lat][:],data[cn_lon][:])
    #plt.plot(x,y,'ob--')
    speeds=[]
    
    
                   
    
    
    tstlon=tstx
    tstlat=tsty
    tstx,tsty=m(tstx,tsty)
    #cmap=plt.cm.hot;
    cmap=arrow_cmap;
    #Q=plt.quiver(tstx,tsty,tstu,tstv,tstc,width=0.0040,headwidth=6,scale=60, color="k",cmap=cmap)
    #THE USUAL SET
    #Q=plt.quiver(tstx,tsty,tstu,tstv,tstc,width=0.004,headwidth=6,scale=60, color="k",cmap=cmap)
    #THE one with 2 cm/s 
    scale_vec=2  #cm/s
    #Q=plt.quiver(tstx,tsty,tstu,tstv,tstc,width=0.004,headwidth=6,scale=60, color="k",cmap=cmap)
    print type(tstx),type(tsty),type(tstu),type(tstv), "TOOOOOTT"
#    print tstx.shape,tsty.shape,tstu.shape,tstv.shape, "TOOOOOTT"
    Q=plt.quiver(tstx,tsty,tstu,tstv,tstc,width=0.004,headwidth=6,scale=15, color="k",cmap=cmap,zorder=1000)
    
    plt.quiverkey(Q,0.9,0.8,scale_vec,str(scale_vec)+r'$ \frac{cm}{s}$',fontproperties={'size':20})
    if (plot_type=='point_number'):
        mn=min_members_in_group
        mx=21
        md=1
        cbar=plt.colorbar(Q,cmap=cmap,boundaries=range(mn,mx,md))
        cbar.set_ticks(range(mx))
    #    cbar.set_ticklabels([mn,md,mx])
        name_base='gridded_argo_surface_amounts'
        if(p_switch=='lower'):
            name_base='gridded_argo_bottom_amounts'
        if(p_gludge_switch=='lower' or p_gludge_switch=='higher'):
            name_base+='_xtra'
    plt.savefig(name_base+'.png' ,facecolor='w',dpi=300)
    plt.savefig(name_base+'.ps' ,facecolor='w',dpi=300)
    plt.savefig(name_base+'.eps' ,facecolor='w',dpi=300)
    print 'average speed for the vectors:', np.average(recorded_speeds)
    print 'Saved figure:', name_base
    
    plt.show()
    vec_data=np.ones((1,5))
    first=True 
    for i in range(len(gridded_vecs)):
        grvec=gridded_vecs[i]
        ucomp=(grvec[2]*lon2m*100)/(grvec[4]) #cm/s
        vcomp=(grvec[3]*lat2m*100)/(grvec[4]) #cm/s
        newline=np.ones((1,5))
        newline[0,:]=np.array([grvec[5],grvec[0],grvec[1],ucomp,vcomp])
        if(first):
            vec_data=np.array(newline)
            first=False
        else:
            vec_data=np.concatenate((vec_data,newline),0)
    np.savetxt(name_base+'.csv',vec_data,'%.2f',',')
    
    #organize the component list
    for j in range(len(orig_vec_data)-1):
        for i in range(j,len(orig_vec_data)-1):
            tmp=list(orig_vec_data[i])
            if(orig_vec_data[i][0]>orig_vec_data[i+1][0]):
                orig_vec_data[i]=list(orig_vec_data[i+1])
                orig_vec_data[i+1]=list(tmp)
            
    #save the component statistic:    
    com_file=open(name_base+'components.csv','w')
    for i in orig_vec_data:
        for comp in i:
            if(type(comp)==int):
                com_file.write("{}".format(comp))
            else:
                if(type(comp)!=str):
                    com_file.write("{:.3f}".format(comp))
                else:
                    com_file.write(comp)
            if(comp!=i[-1]):
                com_file.write(",")
        com_file.write("\n")
            
    com_file.close()
    
argo_to_grid(p_switch='lower',p_gludge_switch='none')
argo_to_grid(p_switch='higher',p_gludge_switch='none')
argo_to_grid(p_switch='lower',p_gludge_switch='higher')
argo_to_grid(p_switch='higher',p_gludge_switch='lower')