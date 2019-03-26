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

max_surface_salinity=7.5
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
    
def is_broken(dat,depth,data_type='salt',qc_flag=None,qc_profs=None):
    dat_change=np.diff(dat)
    d_change=np.diff(depth)
    change=dat_change/d_change
        
    #too high surface salinity
    if(dat[0]>max_surface_salinity):
        return True



    nb_size=6
    #määrittelläänpä pohja oikein.
    true_bottom=len(dat)
    for i in range(len(dat)):
        if not isinstance(dat[i],np.float32):
            true_bottom=i
            break;
    true_bottom-=nb_size  #to avoid checking the very bottom
    for i in range(true_bottom): #don't start at the very bottom, to get number sstabilized
        trigger_jump=0.4
        acceptable_average=0.05
        nb_start=max(0,i-nb_size/2)
        nb_end=min(true_bottom,i+nb_size/2)
#        nb=list(dat_change[nb_start:i])+list(dat_change[i+1:nb_end])
        nb_down=list(dat_change[i+1:nb_end])
        nb_up=list(dat_change[nb_start:i])
        if(len(nb_up)>0 and len(nb_down)>0):
            nb_average_down=sum(nb_down)/len(nb_down)
            nb_average_up=sum(nb_up)/len(nb_up)
            if(np.abs(dat_change[i])>trigger_jump and nb_average_down<acceptable_average):
                print "profile broken at depth:", depth[i],change.shape
                print dat_change[i],nb_average_down, nb_average_up
                return True     #THis signifies a big change in salinitiy in one step,
                            #While rest of the steps hve much smaller steps

    #too big speed of change compared to one of the neighboing points
    """
    accept_ch=0.5  
    for i in range(1,change.shape[0]-1):
        if((np.abs(change[i-1])<np.abs(change[i])*accept_ch or \
            np.abs(change[i+1])<np.abs(change[i])*accept_ch) \
            and np.abs(change[i])>0.2):
            return True     #this means that between two measurements salinity 
                            #has dropped/increased over 0.5 psu, yet
                            #in previous and eralier measurements the change
                            #has been less than tenth of that.
    """

    return False  #at the moment no other checks, thank you.
    
    if(qc_profs is not None):
        for i in qc_profs:
            if(i != '1' and i !=' '):
                print i
                return True
                
        
#    return False

    if(qc_flag is not None):
        if(qc_flag != 'A'):
            print qc_flag
            return True

    
    return False

lineset="2"  #which extra lines to draw. bubblegun/tripwire/ducttape
new_figure=False
value='salt'
#value='oxygen'
current_set=-1
set_no=None
time_min=None
time_max=None

tmin=None
tmax=None
#ref_point = (lats[hsp[6]-1],lons[hsp[6]-1])
ref_point = (57.212,19.806)
ref_dist=100
"""
files=["IM_6902014_20130814_20140821.nc", \
       "IM_6902019_20140821_20150805.nc", \
       "IM_6902020_20150805_20160331_active.nc"]
"""
files=["6902020_20161123123226453.nc",\
       "IM_6902019_20140821_20150805.nc", \
       "IM_6902014_20130814_20140821.nc" \
       
       ]
  
extra_prof_num=0 #add profile number for allready loaded profiles to this. from previous files
prof_no=0
max_prof_no=0
invalid_profiles=0       
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

#figure out how many profiles pass the crieria:
max_prof_no=0
for file_n in files:
    fmk=netcdf.netcdf_file(file_n,'r')
    press=fmk.variables['PRES'][:].copy()
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
    else:
        max_prof_no+=1
    print "max prof no: ", max_prof_no
    fmk.close()




for file_n in files:
    current_set+=1
    fmk=netcdf.netcdf_file(file_n,'r')
    temp=fmk.variables['TEMP'][:].copy()
    salt=fmk.variables['PSAL'][:].copy()
    salt_qc=fmk.variables['PROFILE_PSAL_QC'][:].copy()
    salt_qc_pp=fmk.variables['PSAL_QC'][:].copy()

    press=fmk.variables['PRES'][:].copy()
    oxyg=fmk.variables['DOXY'][:].copy()
    scat=fmk.variables['SCATTERING'][:].copy()

    reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
    jultime = fmk.variables['JULD'][:].tolist()
    apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
    
    
    mask=press>9000
    if(ref_point is not None and ref_dist is not None):
        #We'll limit shown data by distance of the given point
        latits=fmk.variables['LATITUDE'][:].copy()
        longits=fmk.variables['LONGITUDE'][:].copy()
  #      for i in range(np.size(lats)):
  #          if(distance(ref_point,(latits[i],longits[i]))>ref_dist):
  #              mask[i]=True
        
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
    qc_flags=None
    if kuva=='salt_a':
        z_source=salt_a;y_source=press_a;tt=apetime_a;title_txt='PSU'
        qc_flags=salt_qc_a
        qc_points=salt_qc_pp_a
    if kuva=='salt_b':
        z_source=salt_b;y_source=press_b;tt=apetime_b;title_txt='PSU'
        qc_flags=salt_qc_b
        qc_points=salt_qc_pp_b

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
    for i in range(dat.shape[1]-1,-1,-1):
#        if(mask[i,-1] is not True):
         if(distance(ref_point,(latits[i],longits[i]))<=ref_dist):
#            color_profile_d="#{:02x}{:02x}{:02x}".format(min(255,max(0,int(255*(float(prof_no)/max_prof_no)))), \
#                                               min(255,max(0,255-int(255*(float(prof_no)/max_prof_no)))), \
#                                               0*int(255*(float(prof_no)/max_prof_no)))
            c_r=max(min(int(0.25*255*(1.0-float(prof_no)/max_prof_no)),255),0)
            c_g=max(min(int(0.5*255*(1.0-float(prof_no)/max_prof_no)),255),0)
            c_b=max(min(int(255*(1.0-float(prof_no)/max_prof_no)),255),0)
            color_profile_d="#{:02x}{:02x}{:02x}".format(c_r, c_g, c_b)
            alpha=1
            print type(dat),type(pre),type(value),type(qc_flags),type(qc_points)
            if(isinstance(qc_flags,np.ndarray)):
                qc_flag_tmp=qc_flags[i]
            else:
                qc_flag_tmp=None
            
            if(is_broken(dat[:,i],pre[:,i],value,qc_flag_tmp,qc_points[i])):
                invalid_profiles+=1
                color_profile_d="-r"
                alpha=1.0
            else:
                color_profile_d="-k"
                alpha=0.05
            ax.plot(dat[:,i],pre[:,i],color_profile_d,alpha=alpha)
            if value=='oxygen':
                ax2.plot(dat[:,i],pre[:,i],color_profile_d)
            prof_no+=1
    extra_prof_num+=dat.shape[1]
if value=='salt':    
    if(lineset=='1'):
        ax.plot([11,11],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=1)
        ax.plot([12,12],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=1)
        ax.plot([plt.xlim()[0],plt.xlim()[1]],[120,120],color="#808080",linewidth=1)
    if(lineset=='2'):
        ax.plot([max_surface_salinity,max_surface_salinity],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=1)
    plt.ylabel('Pressure [dbar]')
    plt.xlabel('PSU')
if value=='oxygen':    
    plt.axes(ax)
    ax.plot([0,0],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=1)
    ax.plot([15,15],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=1)
    ax.plot([30,30],[plt.ylim()[0],plt.ylim()[1]],color="#808080",linewidth=1)
    ax.plot([plt.xlim()[0],plt.xlim()[1]],[120,120],color="#808080",linewidth=1)
    plt.axes(ax2)
    ax2.plot([plt.xlim()[0],plt.xlim()[1]],[120,120],color="#808080",linewidth=1)
    ax.set_ylim(0,250)
    ax2.set_ylim(0,250)
    ax.set_xlim(0,50)
    ax2.set_xlim(180,400)
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
    plt.xlabel('Oxygen $\mu mol/kg$',labelpad=20)
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
print "total profiles:", prof_no
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