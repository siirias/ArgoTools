# -*- coding: utf-8 -*-
"""
Created on Tue Apr 05 12:39:17 2016

@author: siirias
"""
import sys
sys.path.insert(0,'D:\\svnfmi_merimallit\\qa\\nemo')

import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import ModelQATools as qa
import ModelPltTools
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
runfile('plot_full_data.py')

start=mp.dates.datetime.datetime(1000,5,5)
end=mp.dates.datetime.datetime(3030,5,5)

a1=qa.PointData("6902014_prof.nc",1,start,end,"argonc");
a2=qa.PointData("6902019_prof.nc",1,start,end,"argonc");
a3=qa.PointData("6902020_prof.nc",1,start,end,"argonc");
#a1=qa.PointData("IM_6902014_20130814_20140821.nc",1,start,end,"argonc");
#a2=qa.PointData("IM_6902019_20140821_20150805.nc",1,start,end,"argonc");
#a3=qa.PointData("IM_6902020_20150805_20160331_active.nc",1,start,end,"argonc");


for dataset in range(3):
    if(dataset==0):
        a=a1;offset=0
    if(dataset==1):
        a=a2;offset=a1.obs['ape']['sal'].shape[0]
    if(dataset==2):
        a=a3;offset=a1.obs['ape']['sal'].shape[0]+a2.obs['ape']['sal'].shape[0]
    
    
    data_size=a.obs['ape']['sal'].shape
    print data_size
    print a.obs['ape'].keys()
    d=a.obs['ape']
    #plt.clf();plt.imshow(np.rot90(a.obs['ape']['sal'],3))
    #fig=plt.figure(figsize=(10,10))
    total_im_num=data_size[0]
    total_im_num=1;
    for i in range(total_im_num):
        try:
            plt.close(fig)
        except:
            fir=None
        fig=plt.figure(figsize=(20,10))
        plt.clf()
        
        ax1=fig.add_subplot(121)
    #    ax1.plot(d['lon'][:],-1*d['lat'][:],'k-')
        bmap = Basemap(llcrnrlon=18,llcrnrlat=55,urcrnrlon=21,urcrnrlat=59, resolution = 'i')
        bmap.drawcoastlines()
        bmap.fillcontinents()
        bmap.drawparallels(np.arange(50.,69,1.),labels=[1,0,0,0],linewidth=0)
        bmap.drawmeridians(np.arange(12.,30,1.),labels=[0,0,0,1],linewidth=0)
        x,y=bmap(a1.obs['ape']['lon'][:],a1.obs['ape']['lat'][:])
        bmap.plot(x,y,'k-',linewidth=2)
        x,y=bmap(a2.obs['ape']['lon'][:],a2.obs['ape']['lat'][:])
        bmap.plot(x,y,'g-',linewidth=2)
        x,y=bmap(a3.obs['ape']['lon'][:],a3.obs['ape']['lat'][:])
        bmap.plot(x,y,'b-',linewidth=2)
        #map.plot(np.array([16,21]), np.array([55,59]),'k-',latlon=True)
    #    ax1.plot(a1.obs['ape']['lon'][:],-1*a1.obs['ape']['lat'][:],'k-')
    #    ax1.plot(a2.obs['ape']['lon'][:],-1*a2.obs['ape']['lat'][:],'k-')
    #    ax1.plot(a3.obs['ape']['lon'][:],-1*a3.obs['ape']['lat'][:],'k-')
        
    
        koko=15
    #    plt.plot(d['lon'][i],-1*d['lat'][i],'r*',markersize=koko)
        for tail in range(max(min(koko,i),1)-1,-1,-1):
            x,y=bmap(d['lon'][:],d['lat'][:])
            ax1.plot(x[i-tail],y[i-tail],'ro',markersize=koko-tail,alpha=1.0-float(tail)/koko)
        ax2=fig.add_subplot(222)
        if(i>0):
            for tail in range(max(min(koko,i),1)-1,-1,-1):
                ax2.plot(d['sal'][i-tail,:],d['depth'][i-tail,:],'b--',alpha=1.0-float(tail)/koko)
                ax2.plot(d['tem'][i-tail,:],d['depth'][i-tail,:],'r--',alpha=1.0-float(tail)/koko)

        ax2.plot(d['sal'][i,:],d['depth'][i,:],'b-')
        ax2.plot(d['tem'][i,:],d['depth'][i,:],'r-')

        ax2.set_ylim([0,150])
        ax2.set_xlim([0,20])
        ax2.invert_yaxis()
        plt.ylabel('Pressure [dbar]')
        plt.xlabel('T [$^\circ$C] (Red)/PSU (Blue)')


            
    #    plt.plot(d['sal'][0,:],-1*d['depth'][0,:])
        ax3=plt.subplot(224)
        plot_full_data(set_no=dataset,high_no=i,value="temp")
        """
        plt.set_cmap('Blues')
        ax3.set_axis_bgcolor('black')
        
        var='sal'
        depth=max([a1.obs['ape'][var].shape[1],a2.obs['ape'][var].shape[1],a3.obs['ape'][var].shape[1]])
        dd=a1.obs['ape']
        s=dd[var].shape
        end=s[0]
        im=ax3.imshow(np.rot90(dd[var],3),interpolation="None", aspect='auto',extent=[0,end,depth-s[1],depth])
       
        dd=a2.obs['ape']
        s=dd[var].shape
        im=ax3.imshow(np.rot90(dd[var],3),interpolation="None", aspect='auto',extent=[end,end+s[0],depth-s[1],depth])
        end=end+s[0]
    
        dd=a3.obs['ape']
        s=dd[var].shape
        im=ax3.imshow(np.rot90(dd[var],3),interpolation="None", aspect='auto',extent=[end,end+s[0],depth-s[1],depth])
    
        ax3.plot([i+offset,i+offset],[0,30],'r-',linewidth=2)
        fig.colorbar(im)
        """
        plt.savefig('test_%.3d.png' % (i+offset),facecolor=fig.get_facecolor())
