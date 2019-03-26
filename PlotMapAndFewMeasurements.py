# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 16:11:18 2017

@author: siirias
"""
import sys
sys.path.insert(0,'D:\\svnfmi_merimallit\\qa\\nemo')
import datetime as dt
import calendar
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import ModelQATools as qa
import math
import argohelper as ah

#runfile('D:/ArgoData/plot_full_data.py', wdir='D:/ArgoData')
draw_images=True

#km/day range 
day_in_km=3
how_many_shown=10

        


"""
nc Variable	Variable	Units	Description

metavar1	Cruise		
metavar2	Station		
metavar3	Type		
longitude	Longitude	degrees_east	
latitude	Latitude	degrees_north	
metavar4	Bot. Depth	m	
metavar5	Secchi Depth	m	
date_time	Decimal Gregorian Days of the station	days since 2013-01-01 00:00:00 UTC	Relative Gregorian Days with decimal part
			
var1	PRES	db	
var2	TEMP	deg C	
var3	PSAL	psu	
var4	DOXY	ml/l	
var5	PHOS	umol/l	
var6	TPHS	umol/l	
var7	SLCA	umol/l	
var8	NTRA	umol/l	
var9	NTRI	umol/l	
var10	AMON	umol/l	
var11	NTOT	umol/l	
var12	H2SX	umol/l	
var13	PHPH		
var14	ALKY	meq/l	
var15	CPHL	ug/l	
var16	Year (station date)	
"""
#Gotlands deep
#lon_min=17;lat_min=56;lon_max=22;lat_max=59;
lon_min=16.5;lat_min=55.5;lon_max=22.5;lat_max=59.5;
target_lat=57.3; target_lon=20; target_rad=1.8*6 #rad in km

filetype='csv' # 'nc' tai 'csv'

if filetype=='nc':
    invalid_val=-10000000000.000
    file_n='./Siiriaetal2017/2013-2016_GotlDeep_data_from_helcom.nc'
    fmk=netcdf.netcdf_file(file_n,'r')
    press=-1.0*fmk.variables['var1'][:]
    longitude=fmk.variables['longitude'][:]
    latitude=fmk.variables['latitude'][:]
    start_epoch=dt.datetime(2013,1,1,0,0)
    start_secs=(start_epoch-dt.datetime.utcfromtimestamp(0.0)).total_seconds() #this should be amount of seconds to add to actual timestamps
    #a=dt.datetime.fromtimestamp(times[0]*24.0*60.0*60.0)
    times_s=fmk.variables['date_time'][:]
    times=[]
    for i in range(len(times_s)):
        times.append(dt.datetime.utcfromtimestamp(times_s[i]*24.0*60.0*60.0+start_secs))
    
    d=fmk.variables['var2'][:]
    d=np.ma.masked_where(d==invalid_val,d)
    
else:
    if filetype=='csv':
        invalid_val=-10000000000.000
        #file_n='./Siiriaetal2017/uudet_CTDt_16102017.csv'
        file_n=ah.ctd_data_file
        fmk=pd.read_csv(file_n)
        press=-1.0*fmk[u'PRES [db]'][:]
        longitude=fmk[u'Longitude [degrees_east]'][:]
        latitude=fmk[ u'Latitude [degrees_north]'][:]
        temperature=fmk[u'TEMP [deg C]'][:]
        salinity=fmk[u'PSAL [psu]'][:]
        oxygen=fmk[u'DOXY [ml/l]'][:]
        #purkka, koska happi eri muodossa:
        for i in range(len(oxygen)):
            if(type(oxygen[i])==str and oxygen[i][0]=='<'):
               oxygen[i]=None 
            else:
               oxygen[i]=float(oxygen[i])

        times_s=fmk[u'yyyy-mm-ddThh:mm'][:]
        times=[]
        for i in range(len(times_s)):
            times.append(dt.datetime.strptime(times_s[i],'%Y-%m-%dT%H:%M'))
        ts=[]
        for i in range(len(times)):
            ts.append(calendar.timegm(times[i].timetuple()))
        ctd_data=ah.split_csv_profiles(press,[longitude,latitude,temperature,salinity,ts,oxygen])
       #tehdän se perus aikajana
        timesx=ctd_data[:,0,5]
        times=[]
        for i in range(len(timesx)):
            times.append(dt.datetime.utcfromtimestamp(timesx[i]))  
#        d=np.ma.masked_where(d==invalid_val,d)
        #create mask, for values close to main point:
        distance_mask=[False]*len(times_s)
        for i in range(len(times_s)):
            if(target_rad>ah.distance((target_lat,target_lon),(latitude[i],longitude[i]))):
                distance_mask[i]=True
        #Ja perus paikkajanatkin
        pressure=ctd_data[:,:,0]
        latitude=ctd_data[:,0,2]
        longitude=ctd_data[:,0,1]
        temperature=ctd_data[:,:,3]
        salinity=ctd_data[:,:,4]
        oxygen=ctd_data[:,:,6]
    else:
        print "wrong fieltype,",filetype," aborting!"
        filetype='exitnow'


if(filetype!='exitnow'):
    #
    #
    #
    #
    #
    #create mask, for values close to main point:
    """
    distance_mask=[False]*len(times_s)
    for i in range(len(times_s)):
        if(target_rad>distance((target_lat,target_lon),(latitude[i],longitude[i]))):
            distance_mask[i]=True
    """
    distance_mask=[False]*len(times)
    for i in range(len(times)):
        if(target_rad>ah.distance((target_lat,target_lon),(latitude[i],longitude[i]))):
            distance_mask[i]=True


    
    
    #
    #
    #
    #
    #Map...
    
    if draw_images==True:
        plt.figure(figsize=(15,15))
        bmap = Basemap(projection='merc',\
             resolution='i', llcrnrlat=lat_min, llcrnrlon=lon_min, urcrnrlat=lat_max, urcrnrlon=lon_max)
        bmap.drawcoastlines()
        bmap.fillcontinents()
        bmap.drawparallels(np.arange(50.,69,0.5),labels=[1,0,0,0],linewidth=1,dashes=[1,1],color='#bbbbbb',zorder=0)
        bmap.drawmeridians(np.arange(12.,30,1.0),labels=[0,0,0,1],linewidth=1,dashes=[1,1],color='#bbbbbb',zorder=0)
#        for i in range(len(x)):
#            if(distance_mask[i]):
#                bmap.plot(x[i],y[i],'xg',markersize=15)
        def format_coord(x, y):
            return 'x=%.4f, y=%.4f'%(bmap(x, y, inverse = True))
        plt.gca().format_coord = format_coord
    
    #
    #
    #
    #
    #
    #then plot the argo routes
    plot_legends=False
    plot_routes=True
    
    print "PLOTTING ARGO ROUTES!"
    
    
#    files_to_plot=["6902014_prof.nc","6902019_prof.nc", \
#                        "6902020_prof.nc"]
    files_to_plot=ah.file_names_cleaned
#    colors=["#ff0000","#00ff00","#0000ff"]
    colors=["#5555ff"]*3
    
    labels=[\
            "6902020", "6902021", "6902022", \
            ]
                   
#    start=mp.dates.datetime.datetime(1000,5,5)
#    end=mp.dates.datetime.datetime(3030,5,5)
    start=mp.dates.datetime.datetime(2013,10,19)
    end=mp.dates.datetime.datetime(2013,11,4)
    
    argo_lats=np.array([])
    argo_lons=np.array([])
    argo_times=np.array([])
    argo_depth=[]
    argo_tem=[]
    argo_sal=[]
    
    
    argo_lats_all=np.array([])
    argo_lons_all=np.array([])
    argo_times_all=np.array([])
    argo_depth_all=[]
    argo_tem_all=[]
    argo_sal_all=[]
    
    if plot_routes:
        for f,col,lab in zip(files_to_plot,colors,labels):
            a=qa.PointData(f,1,start,end,"argonc")
            print "File {} has {} profiles".format(f,len(a.obs['ape']['lat'][:]))
            argo_lats=np.concatenate((argo_lats,a.obs['ape']['lat'][:]))
            argo_lons=np.concatenate((argo_lons,a.obs['ape']['lon'][:]))
            argo_times=np.concatenate((argo_times,mp.dates.num2date(a.obs['ape']['date'][:])))
            tmp_val=a.obs['ape']['tem'][:]
            for i in range(tmp_val.shape[0]):
                argo_tem.append(tmp_val[i,:])
            tmp_val=a.obs['ape']['sal'][:]
            for i in range(tmp_val.shape[0]):
                argo_sal.append(tmp_val[i,:])
            tmp_val=a.obs['ape']['depth'][:]
            for i in range(tmp_val.shape[0]):
                argo_depth.append(tmp_val[i,:])
            
            argo_lats_all=np.concatenate((argo_lats_all,argo_lats))
            argo_lons_all=np.concatenate((argo_lons_all,argo_lons))
            argo_times_all=np.concatenate((argo_times_all,argo_times))
            argo_tem_all+=argo_tem
            argo_sal_all+=argo_sal
            argo_depth_all+=argo_depth
            if(draw_images):
                x,y=bmap(a.obs['ape']['lon'][:],a.obs['ape']['lat'][:])
                bmap.plot(x,y,'x-',markersize=4,color=col,linewidth=2, alpha=0.8)
                bmap.plot(x[0],y[0],'o',markersize=5,color=col,linewidth=2, alpha=0.8)
            #bmap.plot(x[-1],y[-1],'o',color=col,markersize=10,alpha=1.0,label=lab)
    
    if(draw_images):
        x,y=bmap(longitude,latitude)
        bmap.plot(x,y,'*',markersize=4,color='#773300')
        #add BY15 place:
        #BY15 (57.32\,N 20.05\,E)
        x,y=bmap([20.05],[57.32])
        bmap.plot(x,y,'o',markersize=10, markerfacecolor='None',color='#ffaaaa')
        plt.savefig('CTD_Argo_places.eps',dpi=300)
        #plot 30 km circle
        angle=np.arange(0,2.0*np.pi,0.01)
        x,y=bmap(20.05+np.sin(angle)*0.5,57.32+np.cos(angle)*0.27)
        bmap.plot(x,y,'r-',color='#aaaaff')
    if(draw_images and False):
        ref_point=(57.315,20.04)
        mbi_dat=0
        tolerance=5.0
        plt.figure(figsize=(12,8))
        for i in range(len(latitude)):
            if(ah.distance(ref_point,(latitude[i],longitude[i]))<50.0):
                plt.plot(oxygen[i],pressure[i])
                fill_data={'lat':latitude[i],'lon':longitude[i],'time':pd.to_datetime(str(times[i])),
                                'sal120':ah.get_closest(pressure[i],salinity[i],-120.0,tolerance)[1],
                                'oxy120':ah.get_closest(pressure[i],oxygen[i],-120.0,tolerance)[1],
                                'sal150':ah.get_closest(pressure[i],salinity[i],-150.0,tolerance)[1],
                                'oxy150':ah.get_closest(pressure[i],oxygen[i],-150.0,tolerance)[1],
                                'sal200':ah.get_closest(pressure[i],salinity[i],-200.0,tolerance)[1],
                                'oxy200':ah.get_closest(pressure[i],oxygen[i],-200.0,tolerance)[1]}
                                                          
                if(type(mbi_dat)==int):
                    mbi_dat=pd.DataFrame(fill_data,index=[0])
                else:
                    idx=len(mbi_dat)
                    mbi_dat=mbi_dat.append(fill_data,ignore_index=idx)
#        mbi_dat=mbi_dat.dropna()
        mbi_dat=mbi_dat.sort_values(by='time')
        plt.figure(figsize=(12,8))
        plt.plot(mbi_dat['time'],mbi_dat['oxy120'],'b*')
        plt.plot(mbi_dat['time'],mbi_dat['oxy150'],'r*')
        plt.plot(mbi_dat['time'],mbi_dat['oxy200'],'g*')
        print mbi_dat.shape
                
#        plt.plot([[i]],[longitude[i]],'x',markersize=4,color=col,linewidth=2, alpha=0.8)
                
#        if(type(mbi_dat)==int):
#            pd.DataFrame()


