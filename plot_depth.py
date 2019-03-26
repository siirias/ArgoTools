# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 10:38:38 2016

@author: siirias
"""

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
import argohelper as ah


def plot_depth( value="temp",target_pressure=120, time_highlight=None, \
                         set_no=None, high_no=None, vmin=None, vmax=None, \
                         new_fig=None, save_file=None, plot_contour=False, \
                         ref_point=None, ref_dist=None,tmin=None,tmax=None, color='k'):

    def distance(origin, destination): 
        lat1, lon1 = origin 
        lat2, lon2 = destination 
        radius = 6371 # km 
        dlat = math.radians(lat2-lat1) 
        dlon = math.radians(lon2-lon1) 
        a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
            * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        d = radius * c 
        return d



#    files=["IM_6902014_20130814_20140821.nc", \
#           "IM_6902019_20140821_20150805.nc", \
#           "IM_6902020_20150805_20160331_active.nc"]

    files=["IM_6902014_20130814_20140821.nc", \
           "IM_6902019_20140821_20150805.nc", \
           "IM_6902020_20150805_20160331_active.nc"]
           
    
#    files=["noora_6902014_20160615160609287test.nc", \
#           "noora_6902019_20160614135812934test.nc", \
#           "noora_6902020_20170828072009212test.nc"]
    files = ah.file_names_cleaned
       
    if(new_fig!=None):
        try:
            fig=plt.figure(figsize=new_fig)   
            plt.clf()
        except:
            fig=plt.figure()   
           
    if(time_highlight!=None):
        time_highlight=mp.dates.date2num(datetime.datetime.strptime(time_highlight,"%Y%m%d"))
    if(value=='temp'):
        kuva='temp_b'
    if(value=='salt'):
        kuva='salt_b'
    if(value=='temp_coarse'):
        kuva='temp_a'
    if(value=='salt_coarse'):
        kuva='salt_a'
    if(value=='oxygen'):
        kuva='oxyg_a'
    if(value=='scatter'):
        kuva='scat_a'
    current_set=-1
    time_min=None
    time_max=None
    for file_n in files:
        current_set+=1
        fmk=netcdf.netcdf_file(file_n,'r')
        temp=fmk.variables['TEMP'][:].copy()
#        temp=fmk.variables['TEMP_ADJUSTED'][:].copy()
        salt=fmk.variables['PSAL'][:].copy()
#        salt=fmk.variables['PSAL_ADJUSTED'][:].copy()
        press=fmk.variables['PRES'][:].copy()
#        press=fmk.variables['PRES_ADJUSTED'][:].copy()
        oxyg=fmk.variables['DOXY'][:].copy()
        scat=fmk.variables['SCATTERING'][:].copy()

        reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
#        reftime = datetime.datetime.strptime(fmk.variables['TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
        jultime = fmk.variables['JULD'][:].tolist()
        apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
        
        
        mask=press>9000
        if(ref_point is not None and ref_dist is not None):
            #We'll limit shown data by distance of the given point
            lats=fmk.variables['LATITUDE'][:].copy()
            lons=fmk.variables['LONGITUDE'][:].copy()
            for i in range(np.size(lats)):
                if(distance(ref_point,(lats[i],lons[i]))>ref_dist):
                    mask[i]=True
            
            
        press_m=press[:].copy()
        press_m[mask]=np.nan
        
        temp_m=temp[:].copy()
        temp_m[mask]=np.nan
        
        salt_m=salt[:].copy()
        salt_m[mask]=np.nan
        
        oxyg_m=oxyg[:].copy()
        oxyg_m[mask]=np.nan
        
        scat_m=scat[:].copy()
        scat_m[mask]=np.nan
        
        apetime_a=apetime[::2].copy()
        apetime_b=apetime[1::2].copy()
        press_a=press_m[::2][:].copy()
        press_b=press_m[1::2][:].copy()
    
        salt_a=salt_m[::2][:].copy()
        salt_b=salt_m[1::2][:].copy()
        temp_a=temp_m[::2][:].copy()
        temp_b=temp_m[1::2][:].copy()
        oxyg_a=oxyg_m[::2][:].copy()
        oxyg_b=oxyg_m[1::2][:].copy()
        scat_a=scat_m[::2][:].copy()
        scat_b=scat_m[1::2][:].copy()
        
        if kuva=='salt_a':
            z_source=salt_a;y_source=press_a;tt=apetime_a;title_txt='PSU'
        if kuva=='salt_b':
            z_source=salt_b;y_source=press_b;tt=apetime_b;title_txt='PSU'
    
        if kuva=='temp_a':
            z_source=temp_a;y_source=press_a;tt=apetime_a;title_txt='Temperature [$^\circ$C]'
        if kuva=='temp_b':
            z_source=temp_b;y_source=press_b;tt=apetime_b;title_txt='Temperature [$^\circ$C]'
    
        if kuva=='oxyg_a':
            z_source=oxyg_a;y_source=press_a;tt=apetime_a;title_txt='Oxygen [$\mu mol/kg$]'
        if kuva=='oxyg_b':
            z_source=oxyg_b;y_source=press_b;tt=apetime_b;title_txt='Oxygen [$\mu mol/kg$]'
        
        if kuva=='scat_a':
            z_source=scat_a;y_source=press_a;tt=apetime_a;title_txt='Scattering [$M^{-1} sr^{-1}$]'
        if kuva=='scat_b':
            z_source=scat_b;y_source=press_b;tt=apetime_b;title_txt='Scattering [$M^{-1} sr^{-1}$]'
        
        
        x_source=np.tile(tt,(z_source.shape[1],1))
        
        x,y=np.mgrid[0:z_source.shape[0],0:z_source.shape[1]]
        x=x.T
        y=-1*y.T
        dat=np.ma.masked_invalid(z_source.T) # z_source.T
        pre=np.ma.masked_invalid(y_source.T)
        time=x_source
        
        if(set_no is not None and high_no is not None):
            if(set_no==current_set):
                time_highlight=time[0][high_no]

        #plt.pcolor(x,y,np.ma.masked_invalid(np.rot90(temp_m)))
    #    plt.figure()
    #    plt.pcolor(x,y,dat)
    #    plt.gca().set_axis_bgcolor('gray')
    #   plt.colorbar()
        
    #    plt.figure()

        #ACTUAL PLOT 
        #plt.pcolor(time,pre,dat,cmap=col_map)

        index_ax=[]
        time_ax=[]
        data_set=[]
        pressure_set=[]
        for i in range(pre.shape[1]):
            index=np.abs(pre[:,i]-target_pressure).argmin()
            diff=np.abs(pre[index,i]-target_pressure)
            index_ax.append(index)
            time_ax.append(time[index,i])
            data_set.append(dat[index,i])
            pressure_set.append(pre[index,i])
            if(diff>3): #this depth doesn't exist
                time_ax[i]=np.nan                
                data_set[i]=np.nan                
        plot_result=plt.plot(np.ma.masked_invalid(time_ax),np.ma.masked_invalid(data_set),color)
#        plt.plot(np.ma.masked_invalid(time_ax),np.ma.masked_invalid(pressure_set),color)
        ax=plt.gca()
#        ax.set_axis_bgcolor('gray')
        plt.ylabel(title_txt)
        plt.xlabel('Time [month-year]')
        #plt.ylim((0,250))
        
        
        #plt.title(title_txt)
#        plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
        plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
        plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%m-%y'))
        #plt.set_cmap('gist_stern')
        if (value=="salt") :
             plt.plot(  [plt.xlim()[0],plt.xlim()[1]], \
                        [11,11],color="#a0a0a0",linewidth=1)
             plt.plot(  [plt.xlim()[0],plt.xlim()[1]], \
                        [12,12],color="#a0a0a0",linewidth=1)

        
        if time_highlight is not None:
            plt.plot([time_highlight,time_highlight],[plt.ylim()[0],plt.ylim()[1]] \
            ,color="#ff0000",linewidth=3)
        if(time_min is None or time.min()<time_min):
                time_min=time.min()
        if(time_max is None or time.max()>time_max):
            time_max=time.max()
    if(tmin is not None):
        time_min=tmin
    if(tmax is not None):
        time_max=tmax
    print time_max,time_min
    print mp.dates.num2date(time_min)
    print mp.dates.num2date(time_max)
    plt.xlim((time_min,time_max))
    locs,labels = plt.xticks()
    plt.setp(labels,rotation=45)
    #plt.savefig('%s.png' % (kuva),dpi=300)

    if(save_file!=None):
        plt.savefig(save_file,dpi=300)
    return plot_result
    
#line1,=plot_depth(new_fig=(15,10),value="oxygen",target_pressure=150,color='b')
#line2,=plot_depth(value="oxygen",target_pressure=120,color='g')
#line1,=plot_depth(new_fig=(15,10),value="salt",target_pressure=150,color='b')
#line2,=plot_depth(value="salt",target_pressure=120,color='g')
#marker_time=mp.dates.date2num(datetime.datetime.strptime("20140821","%Y%m%d"))
#plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)
#marker_time=mp.dates.date2num(datetime.datetime.strptime("20150805","%Y%m%d"))
#plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#d88888",linewidth=2)

#plt.legend([line1,line2],['150 m','120 m'])