# -*- coding: utf-8 -*-
"""
Created on Wed Jul 11 10:56:35 2018

@author: siirias
"""

import argohelper as ah
import re
import datetime as dt
import os
import matplotlib as mp
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap, shiftgrid, cm

alku=\
"""
(Nov 04 2015 08:30:33, 1051593 sec) AirSystem()          Air-bladder inflation by-passed. Barometer: [126cnt, 6.5"Hg].
(Nov 04 2015 08:31:09, 1051629 sec) Telemetry()          Unable to open "7126.1511040830" in append mode.
(Nov 04 2015 08:31:10, 1051630 sec) GpsServices()        GPS almanac is current.
(Nov 04 2015 08:31:10, 1051630 sec) GpsServices()        Initiating GPS fix acquisition.
(Nov 04 2015 08:31:17, 1051637 sec) gga()                $GPGGA,083110,6135.6144,N,02127.8632,E,0,00,,,M,,M,,*59
(Nov 04 2015 08:31:26, 1051646 sec) gga()                $GPGGA,083120,6135.6144,N,02127.8632,E,0,00,,,M,,M,,*5A
(Nov 04 2015 08:31:37, 1051657 sec) gga()                $GPGGA,083130,6135.6144,N,02127.8632,E,0,00,,,M,,M,,*5B
(Nov 04 2015 08:31:46, 1051666 sec) gga()                $GPGGA,083140,6135.6144,N,02127.8632,E,0,00,,,M,,M,,*5C
(Nov 04 2015 08:31:56, 1051676 sec) gga()                $GPGGA,083150,6135.6144,N,02127.8632,E,0,00,,,M,,M,,*5D
(Nov 04 2015 08:32:25, 1051705 sec) gga()                $GPGGA,083200,6135.5311,N,02127.7798,E,1,12,0.7,25.3,M,18.1,M,,*73
(Nov 04 2015 08:32:25, 1051705 sec) GpsServices()        Profile 176 GPS fix obtained in 75 seconds.
(Nov 04 2015 08:32:25, 1051705 sec) GpsServices()                  lon     lat mm/dd/yyyy hhmmss nsat
(Nov 04 2015 08:32:25, 1051705 sec) GpsServices()        Fix:   21.463  61.592 11/04/2015 083200   12


"""

loppu=\
"""
(Nov 04 2015 10:48:15, 1059855 sec) AirSystem()          Air-bladder inflation by-passed. Barometer: [125cnt, 6.2"Hg].
(Nov 04 2015 10:48:47, 1059887 sec) GpsServices()        GPS almanac is current.
(Nov 04 2015 10:48:47, 1059887 sec) GpsServices()        Initiating GPS fix acquisition.
(Nov 04 2015 10:49:25, 1059925 sec) gga()                $GPGGA,104900,6049.6496,N,02334.3331,E,1,11,0.8,70.6,M,16.9,M,,*7B
(Nov 04 2015 10:49:25, 1059925 sec) GpsServices()        Profile 176 GPS fix obtained in 38 seconds.
(Nov 04 2015 10:49:25, 1059925 sec) GpsServices()                  lon     lat mm/dd/yyyy hhmmss nsat
(Nov 04 2015 10:49:26, 1059926 sec) GpsServices()        Fix:   23.572  60.827 11/04/2015 104900   11
(Nov 04 2015 10:49:55, 1059955 sec) gga()                $GPGGA,104930,6049.6480,N,02334.3400,E,1,11,0.8,71.5,M,16.9,M,,*78

"""
plotting=True
if(plotting):
    plt.clf()
    lon_min=12;lat_min=53;lon_max=30;lat_max=66;
    bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
    resolution = 'l',fix_aspect=False)
    bmap.drawcoastlines()

datadirectory="C:\\Data\\ArgoData\\NopeusVertailut\\APE22016\\"
all_files=os.listdir(datadirectory)
files=[]
for f in all_files:
    if(re.match(".*log",f) is not None):
        files.append(f)
for i in range(len(files)-1):
    alku="".join(open(datadirectory+files[i],'r').readlines())
    loppu="".join(open(datadirectory+files[i+1],'r').readlines())
    try:
        t1=dt.datetime.strptime(re.search("\d\d/../.... ......",alku).group(0),"%m/%d/%Y %H%M%S")
        t2=dt.datetime.strptime(re.search("\d\d/../.... ......",loppu).group(0),"%m/%d/%Y %H%M%S")
        
        lat1=float(re.search("Fix:[\s]*(-*\d*\.\d*)\s*(-*\d*\.\d*)",alku).group(2))
        lon1=float(re.search("Fix:[\s]*(-*\d*\.\d*)\s*(-*\d*\.\d*)",alku).group(1))
        
        aika=(t2-t1).total_seconds()
        lat2=float(re.search("Fix:[\s]*(-*\d*\.\d*)\s*(-*\d*\.\d*)",loppu).group(2))
        lon2=float(re.search("Fix:[\s]*(-*\d*\.\d*)\s*(-*\d*\.\d*)",loppu).group(1))
        dist=ah.distance((lat1,lon1),(lat2,lon2))
        
#        print aika/60., 'min'
#        print '{:.2f} km'.format(dist)
#        print '{:.2} m/s'.format(1000.*dist/aika) 
        print("Move: {:04.2f} km in \t{:02.2f} min.\tSpeed {:02.2f} cm/s\t({},{}) to ({},{}) ".format(
                dist, aika/60, 100000.*dist/aika, lat1,lon1,lat2,lon2
                ))       
#        print lat1,lon1,lat2,lon2
        if(plotting):
            x,y=bmap([lon1,lon2],[lat1,lat2])
            bmap.plot(x,y,'-*r')
    except AttributeError:
        print("broken file")