# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 12:04:17 2016

@author: siirias
"""
import math
import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
import cmocean
import argohelper as ah


def plot_full_data( value="temp",col_map='cool', time_highlight=None, \
                         set_no=None, high_no=None, vmin=None, vmax=None, \
                         new_fig=None, save_file=None, plot_contour=False, \
                         ref_point=None, ref_dist=None,tmin=None,tmax=None, \
                         use_converted=False, show_colorbar=True, contour_levels=None,\
                         ymin=0, ymax=250, xlabel=None, ylabel=None,background_color='gray', plot_separators=True):

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


    NONUM_DEF=99999.0
    
           
    files=["noora_6902014_20160615160609287test.nc", \
           "noora_6902019_20160614135812934test.nc", \
           "noora_6902020_20170828072009212test.nc"]
    files=ah.file_names_cleaned
#    files=["6902014_20161123144244280.nc", \
#           "IM_6902019_20140821_20150805.nc", \
#           "6902020_20161123123226453.nc"]
           
           #6902014_20161123144244280.nc
           #6902019_20161123144137259.nc
           #6902020_20161123123226453.nc
           
           #IM_6902014_20130814_20140821.nc
           #IM_6902019_20140821_20150805.nc
           #IM_6902020_20150805_20160331_active.nc
           
    if(new_fig!=None):
        try:
            fig=plt.figure(figsize=new_fig)   
            plt.clf()
        except:
            fig=plt.figure()   
           
    if(time_highlight!=None):
        time_highlight=mp.dates.date2num(datetime.datetime.strptime(time_highlight,"%Y%m%d"))
    kuva=value
    if(value=='temp'):
        kuva='temp_b'
    if(value=='salt'):
        kuva='salt_b'
    if(value=='temp_coarse'):
        kuva='temp_a'
    if(value=='salt_coarse'):
        kuva='salt_a'
    if(value=='oxygen'):
        kuva='oxyg_b'
#    if(value=='scatter'):  #This one requires bit more tinkering below
#        kuva='scat_a'
    current_set=-1
    time_min=None
    time_max=None
    for file_n in files:
        current_set+=1
        fmk=netcdf.netcdf_file(file_n,'r')
        if(use_converted==False):
#            temp=fmk.variables['TEMP'][:].copy()
#            salt=fmk.variables['PSAL'][:].copy()
#            press=fmk.variables['PRES'][:].copy()
#            oxyg=fmk.variables['DOXY'][:].copy()
#            scat=fmk.variables['SCATTERING'][:].copy()

#            temp=fmk.variables['TEMP_ADJUSTED'][:].copy()
#            salt=fmk.variables['PSAL_ADJUSTED'][:].copy()
#            press=fmk.variables['PRES_ADJUSTED'][:].copy()
#            oxyg=fmk.variables['DOXY_ADJUSTED'][:].copy()
#            scat=fmk.variables['SCATTERING_ADJUSTED'][:].copy()
            
            temp=fmk.variables['TEMP_ADJUSTED'][:].copy()
            salt=fmk.variables['PSAL_ADJUSTED'][:].copy()
            press=fmk.variables['PRES_ADJUSTED'][:].copy()
            oxyg=fmk.variables['DOXY'][:].copy()
            scat=fmk.variables['SCATTERING'][:].copy()
    
#            reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
#            jultime = fmk.variables['JULD'][:].tolist()
#            apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
            #####apetime = np.array([datetime.datetime.fromordinal(x) for x in fmk.variables['TIME'][:]])
            if 'TIME' in fmk.variables.keys():            
                apetime = fmk.variables['TIME'][:]
            else:
                reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
                jultime = fmk.variables['JULD'][:].tolist()
                apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
            
        else:
            temp=fmk.variables['var26'][:].copy()
            salt=fmk.variables['var19'][:].copy()
            press=fmk.variables['var23'][:].copy()
            oxyg=fmk.variables['metavar1'][:].copy()
            scat=fmk.variables['var32'][:].copy()
    
            reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
            jultime = fmk.variables['JULD'][:].tolist()
            apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
            
        if(value=='oxygen' and file_n==files[1]):
            kuva='oxyg_a' #Himmeä purkka, koska filet on erilaisia...
        if(value=='oxygen' and file_n!=files[1]):
            kuva='oxyg_b' #Himmeä purkka, koska filet on erilaisia...
        
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
        sa_mean=scat_a[~np.isnan(scat_a)].mean()
        sb_mean=scat_b[~np.isnan(scat_b)].mean()
        scat_gludge=scat_a;
        press_gludge=press_a;
        apetime_gludge=apetime_a
        if(sa_mean>NONUM_DEF*0.9):
            scat_gludge=scat_b
            press_gludge=press_b;
            apetime_gludge=apetime_b
        
        
        
        if kuva=='salt_a':
            z_source=salt_a;y_source=press_a;tt=apetime_a;title_txt='Salinity [g/kg]'
        if kuva=='salt_b':
            z_source=salt_b;y_source=press_b;tt=apetime_b;title_txt='Salinity [g/kg]'
    
        if kuva=='temp_a':
            z_source=temp_a;y_source=press_a;tt=apetime_a;title_txt='Temperature [$^\circ$C]'
        if kuva=='temp_b':
            z_source=temp_b;y_source=press_b;tt=apetime_b;title_txt='Temperature [$^\circ$C]'
    
        
        if kuva=='oxyg_a':
            z_source=oxyg_a;y_source=press_a;tt=apetime_a;title_txt='Oxygen [$ml/l$]'
        if kuva=='oxyg_b':
            z_source=oxyg_b;y_source=press_b;tt=apetime_b;title_txt='Oxygen [$ml/l$]'
        #gludge to switch to another set of oxygen values, if the otehr one is completely empty, but the other is not
        if kuva=='oxyg_a' or kuva=='oxyg_b':
            ox_a=float(np.ma.masked_invalid(oxyg_a).sum())
            ox_b=float(np.ma.masked_invalid(oxyg_b).sum())
            if(ox_a>0 and not (ox_b>0)):
                z_source=oxyg_a;y_source=press_a;tt=apetime_a;title_txt='Oxygen [$ml/l$]'
            if(ox_b>0 and not (ox_a>0)):
                z_source=oxyg_b;y_source=press_b;tt=apetime_b;title_txt='Oxygen [$ml/l$]'        
        
        
        if kuva=='scat_a':
            z_source=scat_a;y_source=press_a;tt=apetime_a;title_txt='Scattering [$M^{-1} sr^{-1}$]'
        if kuva=='scat_b':
            z_source=scat_b;y_source=press_b;tt=apetime_b;title_txt='Scattering [$M^{-1} sr^{-1}$]'
        if kuva=='scatter':
            z_source=scat_gludge;y_source=press_gludge;tt=apetime_gludge;title_txt='Scattering [$M^{-1} sr^{-1}$]'
        
        
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
        try:
            if(plot_contour):
                plt.pcolor( time, \
                            pre,  \
                            dat, \
                            cmap=col_map,vmin=vmin,vmax=vmax)
                try:
                    plt.contourf( time, \
                                pre,  \
                                dat, \
                                cmap=col_map,vmin=vmin,vmax=vmax, \
                                levels=contour_levels)
                except:
                    plt.contourf( time, \
                                pre,  \
                                dat, \
                                cmap=col_map,vmin=vmin,vmax=vmax)
                    
            else:
                plt.pcolor( time, \
                            pre,  \
                            dat, \
                            cmap=col_map,vmin=vmin,vmax=vmax)
        except:
            print "tölölöööö"
            pass
        #plt.plot(time,pre,'*') #Poista tämä rivi, tällä voi tarkastaa missä itse mittpisteet ovat
        
        ax=plt.gca()
        ax.set_facecolor(background_color)
        if(ylabel!=None):
            plt.ylabel(ylabel)
        else:            
            plt.ylabel('Pressure [dbar]')
        if(ylabel!=None):
            plt.xlabel(xlabel)
        else:
            plt.xlabel('Time [month-year]')
            
        plt.ylim((ymin,ymax))
        
        
        #plt.title(title_txt)
#        plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
        #plt.set_cmap('gist_stern')
        if plot_separators and current_set>0:
            plt.plot([apetime[0],apetime[0]],[plt.ylim()[0],plt.ylim()[1]] \
            ,color="#ff0000",linewidth=1)
            
        if time_highlight is not None:
            plt.plot([time_highlight,time_highlight],[plt.ylim()[0],plt.ylim()[1]] \
            ,color="#ff0000",linewidth=1)
        if(time_min is None or time.min()<time_min):
                time_min=time.min()
        if(time_max is None or time.max()>time_max):
            time_max=time.max()
    if(tmin is not None):
        time_min=tmin
    if(tmax is not None):
        time_max=tmax
    plt.xlim((time_min,time_max))
    plt.gca().invert_yaxis()
    if(show_colorbar):
        plt.colorbar(label=title_txt)
    locs,labels = plt.xticks()
    plt.setp(labels,rotation=45)
    plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
    plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%m-%y'))
    #plt.savefig('%s.png' % (kuva),dpi=300)
    if(save_file!=None):
        plt.savefig(save_file+'.png',dpi=300)
        plt.savefig(save_file+'.eps',dpi=300)
#        plt.savefig(save_file+'.jpg',dpi=300)
