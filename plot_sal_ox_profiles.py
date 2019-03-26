# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 20:34:52 2016

@author: siirias
"""
import math
import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
import argohelper as ah
surface_salinity_limit=7.5

def distance(origin, destination): 
    lat1, lon1 = origin 
    lat2, lon2 = destination 
    radius = 6378 # km 
    dlat = math.radians(lat2-lat1) 
    dlon = math.radians(lon2-lon1) 
    a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
    c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
    d = radius * c 
    return d
    
def plot_sal_ox_profiles(value="oxygen"):
    new_figure=True
    high_no=None
    #value='oxygen' or salt
    if(value=="oxygen"):
        lineset="2"  #which extra lines to draw. bubblegun/tripwire/ducttape
    else:
        lineset="1"  #which extra lines to draw. bubblegun/tripwire/ducttape
    current_set=-1
    set_no=None
    time_min=None
    time_max=None
    
    tmin=None
    tmax=None
    #ref_point = (lats[hsp[6]-1],lons[hsp[6]-1])
    ref_point = (57.212,19.806)
    ref_dist=10
    """
    files=["IM_6902014_20130814_20140821.nc", \
           "IM_6902019_20140821_20150805.nc", \
           "IM_6902020_20150805_20160331_active.nc"]
    """
    files=["IM_6902020_20150805_20160331_active.nc",\
           "IM_6902019_20140821_20150805.nc", \
           "IM_6902014_20130814_20140821.nc" \
           
           ]
    files=ah.file_names_cleaned
    #files = ['6902014_20161123144244280_testnew.nc',
    #         '6902019_20161123144137259_testnew.nc',
    #         '6902020_20161123123226453_testnew.nc']
    prof_no=0
    max_prof_no=0
    invalid_profiles=0
    broken_indices={}       
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
    if(value=='scatter'):
        kuva='scat_a'
    if(value=='oxygen'):
    #    ax=fig.add_subplot(121)
    #    ax2=fig.add_subplot(122)
    #    fig,full_ax=plt.subplots(1,1)
        if(new_figure):
            fig=plt.figure(figsize=(10,7),facecolor='white') 
        else:
            plt.clf()
            fig=plt.gcf()
        plt.gca().set_visible(False)
        
        #   full_ax=plt.subplot(gs[0:1])
    
        gs= mp.gridspec.GridSpec(1,2,width_ratios=[3,2])
    
        ax=plt.subplot(gs[0])
        ax2=plt.subplot(gs[1])
    
        full_ax=fig.add_subplot(1,1,1,alpha=0.1,axisbg="#532510")
        full_ax.patch.set_alpha(0.0)
        full_ax.set_visible(False)
        full_ax.set_visible(True)
        full_ax.spines['right'].set_visible(False)
        full_ax.spines['left'].set_visible(False)
        full_ax.spines['top'].set_visible(False)
        full_ax.spines['bottom'].set_visible(False)
        full_ax.set_xticks([])
        full_ax.set_yticks([])
        
    
    #    ax=fig.add_subplot(1,3,1)
    #    ax2=fig.add_subplot(1,3,3)
        plt.subplots_adjust(wspace=0.02)
        ax.spines['right'].set_visible(False)
        ax2.spines['left'].set_visible(False)
        ax.yaxis.tick_left()
        ax2.yaxis.tick_right()
        ax.tick_params(labelright='off')
    
    
    
    else:
        #plt.clf()
        if(new_figure):
            fig=plt.figure(figsize=(10,7))
        else:
            plt.clf()
            fig=plt.gcf()
        ax=plt.gca()
    
    #figure out how many profiles pass teh crieria:
    max_prof_no=0
    min_time=np.nan
    max_time=np.nan
    for file_n in files:
        fmk=netcdf.netcdf_file(file_n,'r')
        reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
        jultime = fmk.variables['JULD'][:].tolist()
        apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
        press=fmk.variables['PRES_ADJUSTED'][:].copy()
        press_m=press[:].copy()
        mask=press>9000
        press_m[mask]=np.nan
        dat=press_m[:][:].copy().T
        if(ref_point is not None and ref_dist is not None):
            #We'll limit shown data by distance of the given point
            latits=fmk.variables['LATITUDE'][:].copy()
            longits=fmk.variables['LONGITUDE'][:].copy()
    
            for i in range(dat.shape[1]-1,-1,-1):
                 if(distance(ref_point,(latits[i],longits[i]))<=ref_dist):
                     max_prof_no+=1
                     if(not np.isnan(min_time) and not np.isnan(max_time)):
                         if(apetime[i]>max_time):
                             max_time=apetime[i]
                         if(apetime[i]<min_time):
                             min_time=apetime[i]
                     else:
                             max_time=apetime[i]
                             min_time=apetime[i]
                        
        else:
            max_prof_no+=1
        print "max prof no: ", max_prof_no
        fmk.close()
    
    
    
    
    for file_n in files:
        current_set+=1
        fmk=netcdf.netcdf_file(file_n,'r')
        temp=fmk.variables['TEMP_ADJUSTED'][:].copy()
        salt=fmk.variables['PSAL_ADJUSTED'][:].copy()
        salt_qc=fmk.variables['PROFILE_PSAL_QC'][:].copy()
    #    salt_qc_pp=fmk.variables['PSAL_QC'][:].copy()
        salt_qc_pp=np.nan
        
        press=fmk.variables['PRES_ADJUSTED'][:].copy()
        oxyg=fmk.variables['DOXY'][:].copy()
        scat=fmk.variables['SCATTERING'][:].copy()
    
        reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
        jultime = fmk.variables['JULD'][:].tolist()
        apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
        
        
        mask=press>9000
    #    if(ref_point is not None and ref_dist is not None):
    #        #We'll limit shown data by distance of the given point
    #        latits=fmk.variables['LATITUDE'][:].copy()
    #        longits=fmk.variables['LONGITUDE'][:].copy()
    #        for i in range(np.size(lats)):
    #            if(distance(ref_point,(latits[i],longits[i]))>ref_dist):
    #                mask[i]=True
            
        fmk.close()        
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
        if(np.isnan(salt_qc_pp)):
            salt_qc_pp_a=np.nan
            salt_qc_pp_b=np.nan
        else:
            salt_qc_pp_a=salt_qc_pp[::2][:].copy()
            salt_qc_pp_b=salt_qc_pp[1::2][:].copy()
        salt_qc_a=salt_qc[::2][:].copy()
        salt_qc_b=salt_qc[1::2][:].copy()
    
        temp_a=temp_m[::2][:].copy()
        temp_b=temp_m[1::2][:].copy()
        oxyg_a=oxyg_m[::2][:].copy()
        oxyg_b=oxyg_m[1::2][:].copy()
        scat_a=scat_m[::2][:].copy()
        scat_b=scat_m[1::2][:].copy()
        if kuva=='salt_a':
            z_source=salt_a;y_source=press_a;tt=apetime_a;title_txt='Salinity [g/kg]'
        if kuva=='salt_b':
            z_source=salt_b;y_source=press_b;tt=apetime_b;title_txt='Salinity [g/kg]'
    
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
        plt.axes(ax)
        broken_indices[file_n]=[]
        for i in range(dat.shape[1]-1,-1,-1):
    #        if(mask[i,-1] is not True):
            if(distance(ref_point,(latits[i],longits[i]))<=ref_dist):
    #            color_profile_d="#{:02x}{:02x}{:02x}".format(min(255,max(0,int(255*(float(prof_no)/max_prof_no)))), \
    #                                               min(255,max(0,255-int(255*(float(prof_no)/max_prof_no)))), \
    #                                               0*int(255*(float(prof_no)/max_prof_no)))
                rate=1.0-(apetime[i]-min_time)/(max_time-min_time)
                c_r=max(min(int(0.25*255*(1.0-rate)),255),0)
                c_g=max(min(int(0.5*255*(1.0-rate)),255),0)
                c_b=max(min(int(255*(1.0-rate)),255),0)
                print rate
                color_profile_d="#{:02x}{:02x}{:02x}".format(c_r, c_g, c_b)
                linew=0.5
                ax.plot(dat[:,i],pre[:,i],color_profile_d,linewidth=linew,zorder=rate)
                if value=='oxygen':
                    ax2.plot(dat[:,i],pre[:,i],color_profile_d,linewidth=linew,zorder=rate)
            prof_no+=1
    if value=='salt':    
        if(lineset=='1'):
#            ax.plot([11,11],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=linew)
#            ax.plot([12,12],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=linew)
            ax.plot([plt.xlim()[0],plt.xlim()[1]],[120,120],color="#808080",linewidth=linew)
        if(lineset=='2'):
            ax.plot([surface_salinity_limit]*2,[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=linew)
        plt.ylabel('Pressure [dbar]')
        plt.xlabel('Salinity [g/kg]')
    if value=='oxygen':    
        plt.axes(ax)
    #    ax.plot([0,0],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=1)
    #    ax.plot([15,15],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=1)
    #    ax.plot([30,30],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=1)
        ax.plot([plt.xlim()[0],plt.xlim()[1]],[120,120],color="#808080",linewidth=1)
        ax.plot([plt.xlim()[0],plt.xlim()[1]],[150,150],color="#808080",linewidth=1)
        plt.axes(ax2)
        ax2.plot([plt.xlim()[0],plt.xlim()[1]],[120,120],color="#808080",linewidth=1)
        ax2.plot([plt.xlim()[0],plt.xlim()[1]],[150,150],color="#808080",linewidth=1)
        ax.set_ylim(0,250)
        ax2.set_ylim(0,250)
        ax.set_xlim(0,1.52)
        ax2.set_xlim(4.5,10)
        ax2.invert_yaxis()
        plt.axes(full_ax)
        br_len=.015
        ax.text(1,0,'/',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontsize='large')
        ax.text(1,1,'/',transform=ax.transAxes,horizontalalignment='center',verticalalignment='center',fontsize='large')
        ax2.text(0,0,'/',transform=ax2.transAxes,horizontalalignment='center',verticalalignment='center',fontsize='large')
        ax2.text(0,1,'/',transform=ax2.transAxes,horizontalalignment='center',verticalalignment='center',fontsize='large')
    #    kwargs=dict(transform=ax.transAxes,color='#000000',clip_on=False,linewidth=1)
    #    ax.plot([1-br_len,1+br_len],[-br_len,br_len],**kwargs)
    #    ax.plot([1-br_len,1+br_len],[1-br_len,1+br_len],**kwargs)
    #    kwargs=dict(transform=ax2.transAxes,color='#000000',clip_on=False,linewidth=1)
    #    ax2.plot([-br_len,br_len],[-br_len,br_len],**kwargs)
    #    ax2.plot([-br_len,br_len],[1-br_len,1+br_len],**kwargs)
    #    ax.plot([xlim[1]+br_len*0.1,xlim[1]-br_len*0.1],[ylim[1]-br_len,ylim [1]+br_len],**kwargs)
    #    xlim=ax2.get_xlim()
    #    ylim=ax2.get_ylim()
    #    ax2.plot([xlim[0]+br_len*0.1,xlim[0]-br_len*0.1],[ylim[0]-br_len,ylim [0]+br_len],color="#000000",linewidth=1,clip_on=False)
    #    ax2.plot([xlim[0]+br_len*0.1,xlim[0]-br_len*0.1],[ylim[1]-br_len,ylim [1]+br_len],color="#000000",linewidth=1,clip_on=False)
    #    ax.plot([0,50],[0,250],color="#00a000",linewidth=12,clip_on=False)
    #    ax.plot([xlim[0],xlim [1]],[ylim[1]-br_len,ylim[1]+br_len],color="#00a000",linewidth=12,clip_on=False)
        plt.ylabel('Pressure [$dbar$]',labelpad=40)
        plt.xlabel('Oxygen $ml/l$',labelpad=20)
        #plt.xlim((-5,80))
    #fig.gca().invert_yaxis()
    ax.invert_yaxis()
    #ax=plt.gca()
    #ax.set_axis_bgcolor('white')
    
    #ax.set_title(title_txt)
    print prof_no
    """    
        try:
            if(plot_contour):
                plt.pcolor( np.ma.masked_invalid(time), \
                            np.ma.masked_invalid(pre),  \
                            np.ma.masked_invalid(dat), \
                            cmap=col_map,vmin=vmin,vmax=vmax)
                plt.contourf( np.ma.masked_invalid(time), \
                            np.ma.masked_invalid(pre),  \
                            np.ma.masked_invalid(dat), \
                            cmap=col_map,vmin=vmin,vmax=vmax)
            else:
                plt.pcolor( np.ma.masked_invalid(time), \
                            np.ma.masked_invalid(pre),  \
                            np.ma.masked_invalid(dat), \
                            cmap=col_map,vmin=vmin,vmax=vmax)
        except:
            pass
        
        ax=plt.gca()
        ax.set_axis_bgcolor('gray')
        plt.ylabel('Pressure [dbar]')
        plt.xlabel('Time [month-year]')
        plt.ylim((0,250))
        
        
        #plt.title(title_txt)
        plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
        plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%m-%y'))
        #plt.set_cmap('gist_stern')
        
        if time_highlight is not None:
            plt.plot([time_highlight,time_highlight],[plt.ylim()[0],plt.ylim()[1]] \
            ,color="#ff0000",linewidth=3)
        if(time_min is None or time.min()<time_min):
                time_min=time.min()
        if(time_max is None or time.max()>time_max):
            time_max=time.max()
        """
    if(tmin is not None):
        time_min=tmin
    if(tmax is not None):
        time_max=tmax
    print time_max,time_min
    
    print "discarded profiles:",invalid_profiles
    for file_n in files:
        print "file {}".format(file_n)
        for i in broken_indices[file_n]:
            print i,
        print ""
    
    plt.savefig("{}_profiles.png".format(value),dpi=300)
    plt.savefig("{}_profiles.eps".format(value),dpi=300)
    
    """
    plt.xlim((time_min,time_max))
    plt.gca().invert_yaxis()
    plt.colorbar(label=title_txt)
    locs,labels = plt.xticks()
    plt.setp(labels,rotation=45)
    #plt.savefig('%s.png' % (kuva),dpi=300)
    
    if(save_file!=None):
        plt.savefig(save_file,dpi=300)
    """    
