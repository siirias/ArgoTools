# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 16:11:18 2017

@author: siirias
"""
import sys
sys.path.insert(0,'D:\\svnfmi_merimallit\\qa\\nemo')
import datetime as dt
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np

#runfile('D:/ArgoData/plot_full_data.py', wdir='D:/ArgoData')


time_min=mp.dates.date2num(dt.datetime.strptime("20150101","%Y%m%d"))
time_max=mp.dates.date2num(dt.datetime.strptime("20151225","%Y%m%d"))
datfile=open('Ostergarnsholm_WS_WD_2013-2016.tsv','r')
dat=datfile.readlines()
datfile.close()
header=dat[0]
dat=dat[1:]
time=[]
speed=[]
direct=[]
for i in range(0,len(dat)):
    dat[i]=dat[i].split('\t')
    dat[i][2]= dt.datetime.strptime(dat[i][2],'%Y-%m-%d %H:%M:%S')
    time.append(dat[i][2])
    dat[i][3]=float(dat[i][3].replace(',','.'))
    speed.append(dat[i][3])
    dat[i][4]=float(dat[i][4].replace(',','.'))
    direct.append(dat[i][4])

#sort the bugger:
order=np.argsort(time)
time= map(lambda i:time[i], order)
speed= map(lambda i:speed[i], order)
direct= map(lambda i:direct[i], order)

fig=plt.figure(figsize=(20,10))
plt.clf()
fig.add_subplot(211)
plt.plot(time,speed,'-')
plt.xlim((time_min,time_max))
plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,1)))
plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%m-%y'))
locs,labels = plt.xticks()
plt.setp(labels,rotation=45)

#plt.plot(time,direct,'-')
fig.add_subplot(212)
plot_full_data(value="oxygen", vmin=0,vmax=50,tmin=time_min,tmax=time_max,show_colorbar=False)
