#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 10:18:27 2017

@author: haaviston
"""

from scipy.io import netcdf
import numpy as np
import datetime
import matplotlib as mp
import gsw
import argohelper as ah
#import argodata

def check_profiles(time,pressure,salinity,QC_prof_salt):
    ok_profiles=range(time.shape[0])
    for i in range(time.shape[0]):    
        ok_profiles[i]=not ah.is_broken(salinity[i][:],pressure[i][:],qc_flag=QC_prof_salt[i],data_type='salt',discard_by_flags=False)
    print time.shape,pressure.shape
    return ok_profiles
#Calculating absolute salinity with TEOS    
def abs_suolaisuus(salt_p,lon,lat):
    a_salt = gsw.SA_from_SP_Baltic(salt_p,lon,lat)
    return np.asarray(a_salt)
    
#Calculating conservative temperature with TEOS    
def conservative_temp(a_salt,temp,press):
    c_temp = gsw.CT_from_t(a_salt,temp,press)
    return np.asarray(c_temp)

def potential_temp(asalt,temp,press):
    pot_temp = gsw.pt_from_t(asalt,temp,press)
    return np.asarray(pot_temp)
    
#Calculating water density with TEOS       
def tiheyslaskin(a_salt,temp,press,lon,lat):
    c_temp = conservative_temp(a_salt,temp,press)
    dens = gsw.rho(a_salt,c_temp,press)
    dens = np.ma.filled(dens,np.nan)
    return dens
    
#This function makes of lon-lat data into a grid    
def lonlatgrid(lons,lats,dat):
    lon = np.empty(np.shape(dat))
    lat = np.empty(np.shape(dat))
    for i in range(lon.shape[0]):
        lon[i,:] = lons[i]
        lat[i,:] = lats[i]
    return lon, lat
    
#This function transforms oxygen data from umol/kg to ml/l
def happimuunnin(oxyg,dens):
    ml_kg_oxyg = oxyg/44.661
    ml_l_oxyg = ml_kg_oxyg*dens/1000    
    return ml_l_oxyg

#This function loads argo-data from a netCDF file
def data_loader(filename,filetype='profile'):
    argo = netcdf.netcdf_file(filename,'r')
#    for v in argo.variables:
#        print v
    if filetype == 'profile':
        temp = argo.variables['TEMP_ADJUSTED'][:].copy()
        salt = argo.variables['PSAL_ADJUSTED'][:].copy()
        press = argo.variables['PRES_ADJUSTED'][:].copy()
    elif filetype == 'traj':
        temp = argo.variables['TEMP_ADJUSTED'][:].copy()
        try:
            argo.variables['PSAL'][:].copy()
        except KeyError:            
            salt = np.full(np.shape(temp),np.nan)
            print "No salinity data :|"
        press = argo.variables['PRES'][:].copy()        
    lats = argo.variables['LATITUDE'][:].copy()
    lons = argo.variables['LONGITUDE'][:].copy()
    reftime_orig=argo.variables['REFERENCE_DATE_TIME'][:]
    reftime = datetime.datetime.strptime(argo.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
    jultime = argo.variables['JULD'][:].tolist()
    time = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
   
    QC_press = argo.variables['PRES_ADJUSTED_QC'][:].copy()
    QC_press[QC_press == ' '] = '0'
    QC_press = QC_press.astype(float)
    QC_press[QC_press == 0] = np.nan
    QC_temp = argo.variables['TEMP_ADJUSTED_QC'][:].copy()
    QC_temp[QC_temp == ' '] = '0'
    QC_temp = QC_temp.astype(float)
    QC_temp[QC_temp == 0] = np.nan
    QC_prof_salt=argo.variables['PROFILE_PSAL_QC'][:].copy()
    QC_prof_salt = QC_prof_salt.astype(str)
    
    QC_salt = argo.variables['PSAL_ADJUSTED_QC'][:].copy()
    QC_salt[QC_salt == ' '] = '0'
    QC_salt = QC_salt.astype(float)
    QC_salt[QC_salt == 0] = np.nan
    QC_date = argo.variables['JULD_QC'][:].copy()
    QC_date[QC_date == ' '] = '0'
    QC_date = QC_date.astype(float)
    QC_date[QC_date == 0] = np.nan
    QC_pos = argo.variables['POSITION_QC'][:].copy()
    QC_pos[QC_pos == ' '] = '0'
    QC_pos = QC_pos.astype(float)
    QC_pos[QC_pos == 0] = np.nan
    
    isBAPE = False
    try:
        oxyg = argo.variables['DOXY'][:].copy()
        QC_oxy = argo.variables['DOXY_QC'][:].copy()
        QC_oxy[QC_oxy == ' '] = '0'
        QC_oxy = QC_oxy.astype(float)
        QC_oxy[QC_oxy == 0] = np.nan
        isBAPE = True
    except KeyError:
        oxyg = np.full(np.shape(temp),np.nan)
        print "No oxygen data :|"
    try:
        scat = argo.variables['SCATTERING'][:].copy()
        QC_scat = argo.variables['SCATTERING_QC'][:].copy()
        QC_scat[QC_scat == ' '] = '0'
        QC_scat = QC_scat.astype(float)
        QC_scat[QC_scat == 0] = np.nan
    except KeyError:   
        scat = np.full(np.shape(temp),np.nan)
        print "No turbidity data :|"
#    !!!!! QC-fägit vaikuttaa kuviin!
#    mask_p = np.logical_or(press>9000,QC_press>=3)
#    mask_t = np.logical_or(temp>9000,QC_temp>=3)
#    mask_s = np.logical_or(salt>9000,QC_salt>=3)
    mask_p = np.where(press>9000)#,QC_press>=3)
    mask_t = np.where(temp>9000)#,QC_temp>=3)
    mask_s = np.where(salt>9000)#,QC_salt>=3)
    
    press_m=press[:].copy()
    press_m[mask_p]=np.nan
    temp_m=temp[:].copy()
    temp_m[mask_t]=np.nan
    salt_m=salt[:].copy()
    salt_m[mask_s]=np.nan
    if isBAPE == True:        
#        mask_o = np.logical_and(oxyg>9000,QC_oxy>=3)
#        mask_s = np.logical_and(scat>9000,QC_scat>=3)
        mask_o = np.where(oxyg>9000)#,QC_oxy>=3)
        mask_s = np.where(scat>9000)#,QC_scat>=3)
        oxyg_m=oxyg[:].copy()
        oxyg_m[mask_o]=np.nan
        scat_m=scat[:].copy()
        scat_m[mask_s]=np.nan
    else:
        oxyg_m = oxyg[:]
        scat_m = scat[:]
        
    mask_pos = np.where(lats>9000)#,QC_pos >= 3)
    lats_m=lats[:].copy()
    lons_m=lons[:].copy()
    lons_m[mask_pos]=np.nan
    lats_m[mask_pos]=np.nan
    
#    mask_time = QC_date==4
    apetime = time[:].copy()
#    apetime[mask_time] = np.nan
#    
    wmo_f = filename[18:25]
    wmo = np.full(np.shape(lats_m),wmo_f)
    wmo =[str(w)[0:7] for w in wmo]
#    
#    if filetype == 'profile':        
#        badprofs = np.array([])
#        for prof in range(np.size(press_m,0)):
#            if np.nanmin(press_m[prof,:]) > 50:
#                badprofs = np.append(badprofs,prof)                
#    
#        press_m = scipy.delete(press_m,badprofs,0)
#        temp_m  = scipy.delete(temp_m,badprofs,0)
#        salt_m  = scipy.delete(salt_m,badprofs,0)
#        oxyg_m  = scipy.delete(oxyg_m,badprofs,0)
#        scat_m  = scipy.delete(scat_m,badprofs,0)
#        apetime = scipy.delete(apetime,badprofs,0)
#        lons_m  = scipy.delete(lons_m,badprofs,0)
#        lats_m  = scipy.delete(lats_m,badprofs,0)
#        wmo     = scipy.delete(wmo,badprofs,0)
#        
#    oxyg_m[oxyg_m==99999.] = np.nan
    
#    scat_m[scat_m>10] = np.nan
        
    return temp_m, salt_m, press_m, oxyg_m, scat_m, apetime, lats_m, lons_m, isBAPE, wmo, QC_prof_salt, reftime_orig, jultime
    
files = ['6902014_20161123144244280',
         '6902019_20161123144137259',
         '6902020_20161123123226453'
         ]
remove_faulty_ones=True
if(remove_faulty_ones):
    file_name_mod='_cleaned.nc'
else:
    file_name_mod='_converted.nc'
    
for f in files:
    temp, salt, press, oxyg, scat, time, lats, lons, isBAPE, wmo, QC_prof_salt, reftime, jultime = data_loader(f+'.nc',filetype='profile')
    
    dens = tiheyslaskin(salt,temp,press,lons,lats)
    ok_profiles=check_profiles(time,press,salt,QC_prof_salt)
    if(not remove_faulty_ones): #Mark everything okay, if not removing the faulty ones
        ok_profiles=len(ok_profiles)*[True]
        
    #These files have two profiles per mesurement,    
    #This loop ensures both are eliminated, if one of them fails the check.
    for i in range(0,len(ok_profiles),2):  
        if((not ok_profiles[i]) or (not ok_profiles[i+1])):
            ok_profiles[i]=False
            ok_profiles[i+1]=False
            
    
    broken_profile_num=np.sum(np.array(ok_profiles)==False) #gives the number of broken profiles
    oxyg = happimuunnin(oxyg,dens)
    print "opening",f+file_name_mod,"for writing"
    net = netcdf.netcdf_file(f+file_name_mod,'w')
    net.history = 'Argo-poijujen happidata ml/l muodossa muiden muuttujien kera'
    
    net.createDimension('oxyprof',np.shape(oxyg)[0]-broken_profile_num)
    net.createDimension('oxydepth',np.shape(oxyg)[1])
    net.createDimension('datetime',np.shape(reftime)[0])

#    net.createDimension('oxyprof',np.shape(oxyg)[0])
#    net.createDimension('oxydepth',np.shape(oxyg)[1])
    presu = net.createVariable('PRES_ADJUSTED','f',('oxyprof','oxydepth'))
    presu[:] = press[ok_profiles][:]
    presu.units = 'dbar'
    
#    net.createDimension('oxyprof',np.shape(oxyg)[0])
#    net.createDimension('oxydepth',np.shape(oxyg)[1])
    sca = net.createVariable('SCATTERING','f',('oxyprof','oxydepth'))
    sca[:] = scat[ok_profiles][:]
    sca.units = 'M-1 sr-1'
    
#    net.createDimension('oxyprof',np.shape(oxyg)[0])
#    net.createDimension('oxydepth',np.shape(oxyg)[1])
    lat = net.createVariable('LATITUDE','f',('oxyprof',))
    lat[:] = lats[ok_profiles][:]
    lat.units = 'deg'
    
#    net.createDimension('oxyprof',np.shape(oxyg)[0])
#    net.createDimension('oxydepth',np.shape(oxyg)[1])
    lon = net.createVariable('LONGITUDE','f',('oxyprof',))
    lon[:] = lons[ok_profiles][:]
    lon.units = 'deg'
    
#    net.createDimension('oxyprof',np.shape(oxyg)[0])
#    net.createDimension('oxydepth',np.shape(oxyg)[1])
    tim = net.createVariable('TIME','f',('oxyprof',))
    tim[:] = time[ok_profiles][:]
    tim.units = 'python datenum'
    
    juld= net.createVariable('JULD','f',('oxyprof',))
    juld[:]=np.array(jultime)[ok_profiles]
    juld.units = "days since 1950-01-01 00:00:00 UTC"
    
    ref_date_time=net.createVariable('REFERENCE_DATE_TIME','c',('datetime',))
    ref_date_time[:]=reftime




    oxygen = net.createVariable('DOXY','f',('oxyprof','oxydepth'))
    oxygen[:] = oxyg[ok_profiles][:]
    oxygen.units = 'ml/l'
    
#    net.createDimension('oxyprof',np.shape(oxyg)[0])
#    net.createDimension('oxydepth',np.shape(oxyg)[1])
    tempe = net.createVariable('TEMP_ADJUSTED','f',('oxyprof','oxydepth'))
    tempe[:] = temp[ok_profiles][:]
    tempe.units = 'deg C'
    
#    net.createDimension('oxyprof',np.shape(oxyg)[0])
#    net.createDimension('oxydepth',np.shape(oxyg)[1])
    salin = net.createVariable('PSAL_ADJUSTED','f',('oxyprof','oxydepth'))
    salin[:] = salt[ok_profiles][:]
#    print "Käännetään PSU:t g/kg:ksi"
    for profile_n in range(salin.shape[0]):
        salin[profile_n,:]=abs_suolaisuus(salin[profile_n,:],lon[profile_n],lat[profile_n])
    salin.units = 'g/kg'
    
    QC_prof_salin= net.createVariable('PROFILE_PSAL_QC','c',('oxyprof',))
    QC_prof_salin[:] = QC_prof_salt[ok_profiles][:]
    QC_prof_salin.units = ''
    
    
    net.sync()
    net.flush()
    net.close()
    
    del net
    del oxygen
    del tempe
    del salin
    del presu
    del lat
    del lon
    del sca
    del tim
    
    print 'lukutesti'
    print f+file_name_mod
    fi = netcdf.netcdf_file(f+file_name_mod, 'r')
    print(fi.history)
    timez = fi.variables['TIME']
    print(timez.units)
    print(timez.shape)
    print(timez[-1])
    print(mp.dates.num2date(timez[-1]))
    fi.close()
    del fi