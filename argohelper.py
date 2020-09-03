# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 10:59:09 2018

@author: siirias
"""
import sys
#sys.path.insert(0,'D:\\svnfmi_merimallit\\qa\\nemo')
import datetime as dt
#import calendar
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import netcdf
#from mpl_toolkits.basemap import Basemap, shiftgrid, cm
#import ModelQATools as qa
import math
import gsw
import cmocean 
import re
ctd_data_file='./Siiriaetal2017/0316544c.csv'

file_names = ['6902014_20161123144244280.nc',
         '6902019_20161123144137259.nc',
         '6902020_20161123123226453.nc']
         
file_names_cleaned=['6902014_20161123144244280_cleaned.nc',
         '6902019_20161123144137259_cleaned.nc',
         '6902020_20161123123226453_cleaned.nc']
file_names_converted=['6902014_20161123144244280_converted.nc',
         '6902019_20161123144137259_converted.nc',
         '6902020_20161123123226453_converted.nc']

def abs_suolaisuus(salt_p,lon,lat):
    a_salt = gsw.SA_from_SP_Baltic(salt_p,lon,lat)
    return np.asarray(a_salt)

def compare_profiles_deprecated(orig_depth,orig_data, comp_depth, comp_data):
    #remove masked values
    if type(orig_depth)==np.ma.core.MaskedArray:
        orig_depth=orig_depth[~orig_depth.mask]
    if type(orig_data)==np.ma.core.MaskedArray:
        orig_data=orig_depth[~orig_depth.mask]
    if type(comp_depth)==np.ma.core.MaskedArray:
        comp_depth=comp_depth[~comp_depth.mask]
    if type(comp_data)==np.ma.core.MaskedArray:
        comp_data=comp_depth[~comp_depth.mask]
    orig_depth=np.abs(np.array(orig_depth))
    comp_depth=np.abs(np.array(comp_depth))
    min_depth=np.nanmax([np.nanmin(orig_depth),np.nanmin(comp_depth)])
    max_depth=np.nanmin([np.nanmax(orig_depth),np.nanmax(comp_depth)])

    total_diff=0
    elements=0
    for i in range(len(orig_depth)):
        if(not np.isnan(orig_depth[i]) and  orig_depth[i]>=min_depth and orig_depth[i]<=max_depth):
            closest=np.nanargmin(np.array(np.abs(comp_depth-orig_depth[i])))
            if(comp_depth[closest]>=min_depth and comp_depth[closest]<=max_depth):
                difference=orig_data[i]-comp_data[closest]
                total_diff+=difference*difference
                elements+=1
    if elements>0:
       total_diff/=elements
    return total_diff

def get_closest(depth,variable,target,tolerance=np.nan):
    if(type(depth)==list):
        depth=np.array(depth)
    if(type(variable)==list):
        variable=np.array(variable)
    depth_dist=abs(depth-target)
    closest=np.nanargmin(depth_dist)
    return (depth[closest],variable[closest])
"""    
    closest_i=0
    distance=abs(depth[0]-target)
    if(np.isnan(variable[0])):
        distance=np.nan
    for i in range(len(depth)):
        if(np.isnan(distance) or distance>abs(depth[i]-target)):
            if(not np.isnan(variable[i])):
                distance=abs(depth[i]-target)
                closest_i=i
    if(np.isnan(tolerance) or distance<tolerance):
        return (depth[closest_i],variable[closest_i])
    else:
        return (np.nan, np.nan)
"""
def compare_profiles(orig_depth,orig_data, comp_depth, comp_data): 
    data,depths,orig_data=difference_profile(orig_depth,orig_data, comp_depth, comp_data)           
    total_diff=0
    elements=0
    for i in range(len(depths)):
        difference=data[i]
        total_diff+=difference*difference
        elements+=1
    if elements>0:
       total_diff/=elements
    return total_diff
            
def difference_profile(orig_depth,orig_data, comp_depth, comp_data):
    #remove masked values
    if type(orig_depth)==np.ma.core.MaskedArray:
        orig_depth=orig_depth[~orig_depth.mask]
    if type(orig_data)==np.ma.core.MaskedArray:
        orig_data=orig_data[~orig_data.mask]
    if type(comp_depth)==np.ma.core.MaskedArray:
        comp_depth=comp_depth[~comp_depth.mask]
    if type(comp_data)==np.ma.core.MaskedArray:
        comp_data=comp_data[~comp_data.mask]
    
    ok_points=[]        
    for i in range(len(comp_depth)):
        if(not np.isnan(comp_data[i])):
            ok_points.append(i)
    comp_data=comp_data[ok_points]
    comp_depth=comp_depth[ok_points]

    ok_points=[]        
    for i in range(len(orig_depth)):
        if(not np.isnan(orig_data[i])):
            ok_points.append(i)
    orig_data=orig_data[ok_points]
    orig_depth=orig_depth[ok_points]
    
    #at this point, check which one is scarser, and use that one as the base:
    flipping=1.0
    if(len(orig_depth)>len(comp_depth)):
        tmp_dep=orig_depth[:].copy()
        tmp_dat=orig_data[:].copy()
        orig_depth=comp_depth[:].copy()
        orig_data=comp_data[:].copy()
        comp_depth=tmp_dep.copy()
        comp_data=tmp_dat.copy()
        flipping=-1.0
    
    orig_depth=np.abs(np.array(orig_depth))
    comp_depth=np.abs(np.array(comp_depth))
    min_depth=np.nanmax([np.nanmin(orig_depth),np.nanmin(comp_depth)])
    max_depth=np.nanmin([np.nanmax(orig_depth),np.nanmax(comp_depth)])
    new_depth=[]
    new_data=[]
    new_orig_data=[]
    for i in range(len(orig_depth)):
        if(not np.isnan(orig_depth[i]) and  orig_depth[i]>=min_depth and orig_depth[i]<=max_depth):
            closest=np.nanargmin(np.array(np.abs(comp_depth-orig_depth[i])))
            if(comp_depth[closest]>=min_depth and comp_depth[closest]<=max_depth):
                difference=orig_data[i]-comp_data[closest]
                new_depth.append(orig_depth[i])
                new_data.append(difference)
                new_orig_data.append(orig_data[i])
    return np.array(new_data)*flipping,np.array(new_depth),np.array(new_orig_data)
    #Flipping is needed, if we swithc the date oterway round, to make sure sarser data is the one used for the dpeth values.


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

def split_csv_profiles(pressure, other_vars, invalid_val=-10000000000.000):
    profiles=0
    depths=0
    max_depths=0
    parameter_num=len(other_vars)+1 #as pressure is one val
    for i in range(1,len(pressure)):
        depths+=1
        if(pressure[i]>pressure[i-1]): #Hypattiin seuraavaan profiiliin
            profiles+=1
            if(max_depths<depths):
                max_depths=depths
            depths=0
    result=np.ones((profiles+1,max_depths+10,parameter_num))*invalid_val
    #Fill the data:
    depth=0
    profile=0
    for i in range(1,len(pressure)):
        if(pressure[i]>pressure[i-1]): #Hypattiin seuraavaan profiiliin
            profile+=1
            depth=0
        if(result[profile][depth][0]!=invalid_val):
            if(abs(result[profile][depth][0]-pressure[i])>0.1):
                print("voi kräpylä", result[profile][depth][0]-pressure[i])
        result[profile][depth][0]=pressure[i]
        for j in range(len(other_vars)):
            result[profile][depth][j+1]=other_vars[j][i]
            
        depth+=1

    result=np.ma.masked_where(result==invalid_val,result)
    return result

def is_broken(dat,depth,data_type='salt',qc_flag=None,qc_profs=None,surface_salinity_limit=7.5,discard_by_flags=True,discard_by_diff=True):
    dat_change=np.diff(dat)
    d_change=np.diff(depth)
    change=dat_change/d_change
    if(data_type!='salt'):
        return False  #doesn't work with anything else yet
# check the QC flag for the profile
    if(discard_by_flags):
        if(qc_flag is not None):
            if(qc_flag != 'A'):
                print(qc_flag)
                return True
#This version of the checks the QC flags for each measurement point
#    if(qc_profs is not None):
#        for i in qc_profs:
#            if(i != '1' and i !=' '):
#                print i
#                return True
    #too high surface salinity
    if(dat[0]>surface_salinity_limit):
        return True
        
    #too big speed of change compared to one of the neighboing points
    #return False
        
    if(discard_by_diff):
        accept_ch=0.1  
        for i in range(1,change.shape[0]-1):
            if((np.abs(change[i-1])<np.abs(change[i])*accept_ch or \
                np.abs(change[i+1])<np.abs(change[i])*accept_ch) \
                and np.abs(change[i])>0.3):
                return True     #this means that between two measurements salinity 
                                #has dropped/increased over 0.5 psu, yet
                                #in previous and eralier measurements the change
                                #has been less than tenth of that.
    
    return False

def give_statistics(files_to_use = None):
    #bit deprecated, check rather gather_statistics()
    start=mp.dates.datetime.datetime(1000,5,5)
    end=mp.dates.datetime.datetime(3030,5,5)
    if(not files_to_use):
        files_to_use = file_names_converted
    elif type(files_to_use) == str:
        files_to_use = [files_to_use]
    for filename in files_to_use:
        print(filename,":")
        argo = netcdf.netcdf_file(filename,'r')
        temp = argo.variables['TEMP_ADJUSTED'][:].copy()
        salt = argo.variables['PSAL_ADJUSTED'][:].copy()
        press = argo.variables['PRES_ADJUSTED'][:].copy()
        lats = argo.variables['LATITUDE'][:].copy()
        lons = argo.variables['LONGITUDE'][:].copy()
        reftime_orig=argo.variables['REFERENCE_DATE_TIME'][:]
        reftime = dt.datetime.strptime(argo.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
        jultime = argo.variables['JULD'][:].tolist()
        time = np.array([mp.dates.date2num(reftime + dt.timedelta(days=x)) for x in jultime])
        
        avg_pres=0.0
        max_pres=0.0
        min_pres=10000000.0
        for i in range(press.shape[0]):
            avg_pres+=max(press[i,:])
            if(max_pres<max(press[i,:])):
                max_pres=max(press[i,:])
            if(min_pres>max(press[i,:])):
                min_pres=max(press[i,:])
        avg_pres/=float(press.shape[0])
        print("deployed: {} and recovered: {}".format(str(mp.dates.num2date(time[0])),str(mp.dates.num2date(time[-1]))))
        print("mission time: {}".format(mp.dates.num2date(time[-1])-mp.dates.num2date(time[0])))
        print("profiles: {}".format(temp.shape[0]))
        print("depth avg:{}, min:{},max:{}".format(avg_pres,min_pres,max_pres))
        print("-----\n")
        
def get_primary_indices(dataset):
    #dataset is netcdffile loaded with xarray .opendataset
    sampling_schemes =  dataset['VERTICAL_SAMPLING_SCHEME']
    #See http://www.odip.org/documents/odip/downloads/20/argo-dm-user-manual.pdf
    # page 18, this determines the type of profile, only one primary
    #per profile
    sampling_schemes = map(str,sampling_schemes) #convert to strings
    primaries = list(map(lambda x:'Primary' in x, sampling_schemes))
    #true where primary.
    primaries=np.array(primaries)
    return primaries

def interpolate_data_to_depths(variable, depths, new_depth_axis):
    #assumes data is format (profile_n, level_n)
    d_shape = variable.shape
    new_data = np.zeros((d_shape[0],new_depth_axis.shape[0]))
    for i in range(new_data.shape[0]):
        min_depth = np.nanmin(depths[i,:])
        max_depth = np.nanmax(depths[i,:])
        for d in range(new_data.shape[1]):
            val = get_closest(depths[i,:],variable[i,:],new_depth_axis[d])[1]
            if(new_depth_axis[d]>max_depth or \
               new_depth_axis[d]<min_depth):
                val = np.nan
            new_data[i,d] = val
    return  new_data

def axes_label_from_variable_name(var_name, give_colormap = False):
    #picks typical axis label based on usual parameters plotted
    axes_label = var_name
    colormap = cmocean.cm.tarn
    if var_name in ['TEMP']:
        axes_label = 'Temperature/ $^\circ$C'
        colormap = cmocean.cm.thermal
    if var_name in ['PSAL']:
        axes_label = 'Salinity'
        colormap = cmocean.cm.haline
    if var_name in ['PRES']:
        axes_label = 'Pressure/dbar'
        colormap = cmocean.cm.deep
    if(give_colormap):
        return (axes_label, colormap)
    else:
        return axes_label

def gather_statistics(dataset, filter_bool = slice(None)):
    stats = {}
    time = dataset['JULD'][filter_bool]
    depths = dataset['PRES'][filter_bool]
    depths = list(map(lambda x: np.max(x), depths))
    stats['deployment_lat'] = float(dataset['LATITUDE'][0])
    stats['deployment_lon'] = float(dataset['LONGITUDE'][0])
    distances_total = []
    distances_last = []
    prev_lat = stats['deployment_lat']
    prev_lon = stats['deployment_lon']
    for lat,lon in zip(dataset['LATITUDE'][filter_bool], \
                       dataset['LONGITUDE'][filter_bool]):
        distances_total.append(\
            distance([stats['deployment_lat'],stats['deployment_lon']],\
                     [lat,lon]))
        distances_last.append(\
            distance([prev_lat,prev_lon],\
                     [lat,lon]))
        prev_lat = lat
        prev_lon = lon
    stats['distance_from_origin'] = distances_total
    stats['distance_since_last'] = distances_last
    stats['wmo'] = str(int(dataset['PLATFORM_NUMBER'][0].data))
    
    stats['serial'] = \
        re.search("'(.*)'",\
        str(dataset['FLOAT_SERIAL_NO'][0].data)).groups()[0].strip()
    stats['type'] = \
        re.search("'(.*)'",\
        str(dataset['PLATFORM_TYPE'][0].data)).groups()[0].strip()
        #This is an awful gludge, but couldn't get the string out outherways...
    
    #even more horrible gludge, to get the sensors:
    tmp_main = map(str,dataset['PARAMETER'][0,0,:].data)
    tmp = list(tmp_main)+list(map(str,dataset['PARAMETER'][1,0,:].data))
    tmp = list(map(lambda x: re.search("'(.*)'",x).groups()[0].strip(),tmp))
    tmp = [x for x in tmp if len(x)>0]
    stats['sensors'] = "-".join(set(tmp))  # set removes dublicates

    stats['time_deployed']=pd.to_datetime(time[0].data)
    stats['time_last_profile']=pd.to_datetime(time[-1].data)
    stats['depth_avg']=np.mean(depths)
    stats['depth_min']=np.min(depths)
    stats['depth_max']=np.max(depths)
    times_between = np.array(list(\
                    map(lambda x: float(x),np.diff(time))))/\
                    (1000000000.0*60.0*60.0) # from ns to hours
    stats['times_between']=times_between
    
    #Bit of gludge, but some hardcoded areas:
    stats['area'] = "{}-{}".format(stats['deployment_lat'],\
                                     stats['deployment_lon'])
    if(stats['deployment_lat']>65.0):
        stats['area'] = "Barents Sea"
    elif(stats['deployment_lat']>63.0):
        stats['area'] = "Bay of Bothnia"
    elif(stats['deployment_lat']>60.0):
        stats['area'] = "Bothnian Sea"
    elif(stats['deployment_lat']>58.0):
        stats['area'] = "N.Baltic Proper"
    elif(stats['deployment_lat']>56.3):
        stats['area'] = "Baltic Proper"
    elif(stats['deployment_lon']>17.3):
        stats['area'] = "Gdansk Basin"
    elif(stats['deployment_lon']<17.3):
        stats['area'] = "Bornholm Basin"
    
    #Another gludge to attach nicknames for floats
    nicknames = {\
                 '023-3119':'EAR-2',\
                 'AC0300-19FI001':'',\
                 'AI2600-18FI001':'Arvo1',\
                 'AI2600-19FI001':'EAR-1',\
                 '023-3119':'',\
                 '6710':'BAPE2',\
                 '6711':'BAPE1',\
                 '7191':'BAPE3',\
                 '7958':'HAPE1',\
                 '7959':'PAPE1',\
                 '5088':'APE1',\
                 '5396':'APE2',\
                 '5397':'APE1',\
                 '8540':'BAPE3',\
                 '8541':'HAPE2',\
                 '8543':'PAPE3',\
                 '8348':'CAPE1',\
                 '9568':'BAPE3',\
                 }
    stats['nickname'] = '-'
    if stats['serial'] in nicknames.keys():
        stats['nickname'] = nicknames[stats['serial']]
        
    
    return stats