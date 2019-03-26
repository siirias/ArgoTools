# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:11:08 2016

@author: siirias
"""
import sys
sys.path.insert(0,'D:\\svnfmi_merimallit\\qa\\nemo')

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import ModelQATools as qa
import ModelPltTools
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from netCDF4 import Dataset

#Full baltic
lon_min=16;lat_min=53;lon_max=30;lat_max=66;
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;
#Bothnian Sea
#lon_min=17;lat_min=60;lon_max=22;lat_max=63;

plot_bathymetry=True
plot_legends=True
plot_routes=True

fig=plt.figure(figsize=(8,10))

#fig=plt.figure(figsize=(10,13))
#fig=plt.figure(figsize=(10,10))


plt.clf()
bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = 'f',fix_aspect=False)

bmap.drawcoastlines()
bmap.fillcontinents()
bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])
"""
filest_to_plot=["6901901_prof.nc","6902017_prof.nc","6902013_prof.nc", \
                "6902018_prof.nc","6902014_prof.nc","6902019_prof.nc", \
                "6902020_prof.nc", "6902021_prof.nc", "6902022_prof.nc", \
                "6902023_prof.nc", "6902024_prof.nc", "6902025_prof.nc",\
                "6902026_prof.nc", "6902027_prof.nc", "6902028_prof.nc", "6902029_prof.nc"\
                ]

                
colors=["#ff0000","#00ff00","#0000ff","#000000","#d00000","#d000d0","#005050", \
        "#888888", "#ffff00", "#ff00ff", "#00ffff", "#ffffff", \
        "#aa0055", "#50ff00", "#3060ff", "#0055ff", "#ff22aa" \
        ]

labels=["6901901","6902017","6902013", \
        "6902018","6902014","6902019", \
        "6902020", "6902021", "6902022", "6902023", "6902024",\
        "6902025", "6902026", "6902027", "6902028", "6902029"\
        ]
"""


##2016-2017
#filest_to_plot=[\
                #"6902020_prof.nc", "6902021_prof.nc", "6902022_prof.nc", \
                #"6902023_prof.nc", "6902024_prof.nc", "6902025_prof.nc",\
                #"6902026_prof.nc", "6902027_prof.nc", "6902028_prof.nc",\
                #"6902029_prof.nc"\
                #]
#
                #
#colors=["#ff0000","#00ff00","#0000ff",\
        #"#000000","#d00000","#d000d0",\
        #"#005050", "#888888", "#ffff00",\
        #"#ff00ff"\
        #]
#
#labels=[\
        #"6902020", "6902021", "6902022", \
        #"6902023", "6902024", "6902025", \
        #"6902026", "6902027", "6902028", \
        #"6902029"\
        #]

"""
#2016
filest_to_plot=[\
                "6902020_prof.nc", "6902021_prof.nc", "6902022_prof.nc", \
                "6902023_prof.nc", "6902024_prof.nc"
                ]
                
colors=["#ff0000","#00ff00","#0000ff",\
        "#000000","#d00000" \
        ]

labels=[\
        "6902020", "6902021", "6902022", \
        "6902023", "6902024" \
        ]
"""

filest_to_plot=["6901901_prof.nc","6902013_prof.nc","6902017_prof.nc", \
                "6902018_prof.nc","6902021_prof.nc", "6902022_prof.nc", \
                "6902023_prof.nc","6902025_prof.nc", "6902028_prof.nc", \
                "6902029_prof.nc", "6902014_prof.nc", "6902019_prof.nc", \
                "6902020_prof.nc"]

                
colors=["#ff0000","#00ff00","#0000ff",\
        "#000000","#d00000","#d000d0",\
        "#005050","#888888", "#ffff00",\
        "#ff00ff", "#00ffff", "#ffffff", \
        "#aa0055"]

labels=["6901901","6902013","6902017", \
        "6902018","6902021","6902022", \
        "6902023", "6902025", "6902028",\
        "6902029", "6902014", "6902019",\
        "6902020" ]


#filest_to_plot=["6901901_prof.nc","6902017_prof.nc","6902013_prof.nc", \
#"6902018_prof.nc", "6902021_prof.nc", "6902022_prof.nc"]

#colors=colors[:len(filest_to_plot)]
#labels=["6901901","6902017","6902013", \
#"6902018", "6902021", "6902022"]
               
start=mp.dates.datetime.datetime(1000,5,5)
end=mp.dates.datetime.datetime(3030,5,5)

#TOPOGRAPHY EXPERIMENT
if plot_bathymetry:
    topodata = Dataset('iowtopo2_rev03.nc')
    
    topoin = topodata.variables['Z_WATER'][:]
    lons = topodata.variables['XT_I'][:]
    lats = topodata.variables['YT_J'][:]
    x=np.tile(lons,(lats.shape[0],1))
    y=np.tile(lats,(lons.shape[0],1)).T
    bmap.pcolor(x,y,-1*topoin,cmap='bone_r',vmin=0,vmax=300)
    cb=plt.colorbar()
    cb.ax.invert_yaxis()

if plot_routes:
    for f,col,lab in zip(filest_to_plot,colors,labels):
        a=qa.PointData(f,1,start,end,"argonc");
        lon_dat=a.obs['ape']['lon'][~a.obs['ape']['lon'].mask]
        lat_dat=a.obs['ape']['lat'][~a.obs['ape']['lat'].mask]
        x,y=bmap(lon_dat,lat_dat)
    #    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
        bmap.plot(x,y,color=col,linewidth=2, alpha=0.8)
        bmap.plot(x[-1],y[-1],'x',color=col,markersize=8,alpha=1.0)
        bmap.plot(x[0],y[0],'o',color=col,markersize=6,alpha=1.0,label=lab)
    #    print lab, mp.dates.num2date(a.obs['ape']['date'][0]).date() \
    #             , mp.dates.num2date(a.obs['ape']['date'][-1]).date()

if plot_legends:
#    plt.legend(bbox_to_anchor=(1,0.4),numpoints=1)
    #plt.legend(bbox_to_anchor=(1.0,0.5),numpoints=1)
    plt.legend(loc='lower right',numpoints=1)
plt.savefig('ArgoRoutesArticles.png' ,facecolor='w',dpi=300)
plt.savefig('ArgoRoutesArticles.eps' ,facecolor='w',dpi=300)
