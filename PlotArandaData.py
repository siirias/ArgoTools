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
    if(draw_images):
        for i in range(len(press[:,0])):
            plt.plot(d[i,:]+0.5*i,press[i,:])
    
else:
    if filetype=='csv':
        invalid_val=-10000000000.000
        file_n='./Siiriaetal2017/Uudet_CTDt_16102017.csv'
        fmk=pd.read_csv(file_n)
        press=-1.0*fmk[u'PRES [db]'][:]
        longitude=fmk[u'Longitude [degrees_east]'][:]
        latitude=fmk[ u'Latitude [degrees_north]'][:]
        temperature=fmk[u'TEMP [deg C]'][:]
        salinity=fmk[u'PSAL [psu]'][:]

        times_s=fmk[u'yyyy-mm-ddThh:mm'][:]
        times=[]
        for i in range(len(times_s)):
            times.append(dt.datetime.strptime(times_s[i],'%Y-%m-%dT%H:%M'))
        ts=[]
        for i in range(len(times)):
            ts.append(calendar.timegm(times[i].timetuple()))
        ctd_data=ah.split_csv_profiles(press,[longitude,latitude,temperature,salinity,ts])
        
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
        if(draw_images):
            for i in range(ctd_data.shape[0]):
                plt.plot(temperature[i,:]+0.5*i,pressure[i,:])
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


    #Let's plot times and measurements.
    if draw_images==True:
        plt.figure()
        plt.plot(times,latitude,'*')
    
    
    
    
    #
    #
    #
    #
    #Map...
    
    if draw_images==True:
        plt.figure(figsize=(6,4))
        bmap = Basemap(projection='merc',\
             resolution='i', llcrnrlat=lat_min, llcrnrlon=lon_min, urcrnrlat=lat_max, urcrnrlon=lon_max)
        bmap.drawcoastlines()
        bmap.fillcontinents()
        bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
        bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])
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
    
    
    files_to_plot=["6902014_prof.nc","6902019_prof.nc", \
                        "6902020_prof.nc"]
    files_to_plot=ah.file_names_cleaned
#    colors=["#ff0000","#00ff00","#0000ff"]
    colors=["#5555ff"]*3
    
    labels=[\
            "6902020", "6902021", "6902022", \
            ]
                   
    start=mp.dates.datetime.datetime(1000,5,5)
    end=mp.dates.datetime.datetime(3030,5,5)
    
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
                bmap.plot(x,y,'x',markersize=4,color=col,linewidth=2, alpha=0.8)
            #bmap.plot(x[-1],y[-1],'o',color=col,markersize=10,alpha=1.0,label=lab)
    
    if(draw_images):
        x,y=bmap(longitude,latitude)
        bmap.plot(x,y,'*',markersize=4,color='#773300')
        plt.savefig('CTD_Argo_places.eps',dpi=300)


    if plot_legends and draw_images:
        plt.legend(bbox_to_anchor=(1.0,0.5),numpoints=1)
        
    #make the distance mask for argo floats
    argo_distance_mask=[False]*len(argo_lats)
    for i in range(len(argo_lats)):
        if(target_rad>ah.distance((target_lat,target_lon),(argo_lats[i],argo_lons[i]))):
            argo_distance_mask[i]=True
        
        
    print "CTD measurements\t:",sum(distance_mask)
    print "Argo measurements\t:",sum(argo_distance_mask)
    
    
    #Plot timeline of picked measurements:
    if draw_images==True:
        plt.figure()
        for i in range(len(distance_mask)):
            if(i):
                dist=ah.distance((target_lat,target_lon),(latitude[i],longitude[i]))
                plt.plot(times[i],dist,'*r')
        for i in range(len(argo_distance_mask)):
            if(i):
                dist=ah.distance((target_lat,target_lon),(argo_lats[i],argo_lons[i]))
                plt.plot(argo_times[i],dist,'og')
    
    #Lets find out the closest Argo measurement, for each CTD measurement:
    best_hits=[]
    for CTD in range(len(times)):
        closest_index=0
        closest_dist=-1.
        closest_time=-1.
        closest_dist_km=-1.
        for Argo in range(len(argo_lats_all)):
            dist_km=ah.distance((argo_lats_all[Argo],argo_lons_all[Argo]),(latitude[CTD],longitude[CTD]))
            dist_days= abs(mp.dates.date2num(argo_times_all[Argo])-mp.dates.date2num(times[CTD]))
            dist=dist_km+day_in_km*dist_days
            if(closest_dist<0 or closest_dist>dist):
                closest_dist=dist
                closest_dist_km=dist_km
                closest_time=dist_days
                closest_index=Argo
        best_hits.append([0,0,0,0,0])
        best_hits[CTD][0]=closest_index
        best_hits[CTD][1]=closest_dist
        best_hits[CTD][2]=closest_time
        best_hits[CTD][3]=closest_dist_km
        best_hits[CTD][4]=CTD
#        print "for CTD",CTD,"Argo",best_hits[CTD][0],"\t\td(km){:.2f} (d){:.0f}".format(best_hits[CTD][2],best_hits[CTD][3])
        
#sort based on best value:
for i in range(len(best_hits)):
    for j in range(i,len(best_hits)):
        tmp=best_hits[i]
        if(best_hits[i][1]>best_hits[j][1]):
            best_hits[i]=best_hits[j]
            best_hits[j]=tmp            
avg_km_diff=0.
avg_t_diff=0.
for i in range(how_many_shown):
    print "for CTD",best_hits[i][4],"Argo",best_hits[i][0],"\t\td(km){:.2f} (d){:.1f} (tot){:.1f}".format(best_hits[i][3],best_hits[i][2],best_hits[i][1])
    avg_km_diff+=best_hits[i][3]
    avg_t_diff+=best_hits[i][2]
avg_km_diff/=how_many_shown
avg_t_diff/=how_many_shown


plt.figure()
graph_step=10.0
total_fitness=0.0
for No in range(how_many_shown):
    iargo=best_hits[No][0]
    ictd=best_hits[No][4]
    plt.plot(temperature[ictd,:]+graph_step*No,pressure[ictd,:],'r-')
#    plt.plot(ctd_data[ictd,:,3]+0.5*No,ctd_data[i,:,0],'-x')
    plt.plot(argo_tem_all[iargo][:]+graph_step*No,-1.0*argo_depth_all[iargo][:],'k-')
    fitness=ah.compare_profiles(pressure[ictd,:],temperature[ictd,:],argo_depth_all[iargo][:],argo_tem_all[iargo][:])
    total_fitness+=fitness
    print "profile {} fitness {}".format(No,fitness)
    print "Distance in km: {} \t Distance in days: {}".format(best_hits[No][3],best_hits[No][2])
    print "Location (CTD) {},{}".format(latitude[best_hits[No][4]],longitude[best_hits[No][4]])
    print "Location (Argo) {},{}".format(argo_lats_all[best_hits[No][0]],argo_lons_all[best_hits[No][0]])
    print "Time (CTD) {}".format(times[best_hits[No][4]])
    print "Time (Argo) {}\n\n".format(argo_times_all[best_hits[No][0]])
print "With magic number {}, \taverage fitness: {:.2f} \tavg km diff {:.1f} \tavg d diff {:.2f}".format(day_in_km,total_fitness/how_many_shown,avg_km_diff, avg_t_diff)
print "FINISHED!"