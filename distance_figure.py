# -*- coding: utf-8 -*-
"""
Created on Mon Jun 06 15:38:26 2016

@author: siirias
"""
import numpy as np
import sys
import matplotlib.pyplot as plt
hsp=[33,45,67,98,116,134,186]
#print lats[hsp],lons[hsp]
import plot_full_distance as pfdist


"""
def draw_map(lats,lons,hsp):
        try:
            plt.close(fig)
        except:
            fir=None
        fig=plt.figure(figsize=(8,12))
        plt.clf()    
        bmap = Basemap(llcrnrlon=18,llcrnrlat=55,urcrnrlon=21,urcrnrlat=59, resolution = 'i')
        bmap.drawcoastlines()
        bmap.fillcontinents()
        bmap.drawparallels(np.arange(50.,69,1.),labels=[1,0,0,0],linewidth=0)
        bmap.drawmeridians(np.arange(12.,30,1.),labels=[0,0,0,1],linewidth=0)
        x,y=bmap(a1.obs['ape']['lon'][:],a1.obs['ape']['lat'][:])
        bmap.plot(x,y,'k-',linewidth=2)
        x,y=bmap(a2.obs['ape']['lon'][:],a2.obs['ape']['lat'][:])
        bmap.plot(x,y,'k-',linewidth=2)
        x,y=bmap(a3.obs['ape']['lon'][:],a3.obs['ape']['lat'][:])
        bmap.plot(x,y,'k-',linewidth=2)
        for i in hsp:
            x,y=bmap(lons[i],lats[i])
            plt.plot(x,y,'ro',markersize=15,alpha=0.8)
           #plt.text(x-0.3,y+0.0,'{}'.format(i),color='blue')
"""        
#sys.exit("")
print "Calculating distances"
avg_distances=[]
points=201
max_dist=100
dist_under=np.zeros((points,max_dist))
for i in range(points):
    (time,dist,lats,lons)=pfdist.plot_full_distance(ref_is_no=i,vmin=0,vmax=160,no_plot=True)
    avg_distances.append(np.mean(dist))
    for d in range(max_dist):
        for j in range(points):
            if(d>dist[j]):
                dist_under[i,d]+=1
#    plt.plot([time[i],time[i]],[plt.ylim()[0],plt.ylim()[1]] \
#    ,color="#ff0000",linewidth=3)
        
#    save_file="distances_{:03}.png".format(i)
#    plt.savefig(save_file,dpi=120)
    if(dist_under[i,20]>20):
        print lats[i],lons[i],dist_under[i,20]
"""
print "Drawing nearby points"
target_distance=20
for i in range(points):
    DT=mp.dates.num2date(time[i])
    date_str="{:04}{:02}{:02}".format(DT.year,DT.month,DT.day)
    plt.close();plt.clf()
    plot_full_data(new_fig=True,ref_point=(lats[i],lons[i]), ref_dist=target_distance,time_highlight=date_str ,value="salt", vmin=11,vmax=13,tmin=735600,tmax=735850)
    save_file="Salt_cut_dist_{:02}_{:03}.png".format(target_distance,i)
    plt.savefig(save_file,dpi=120)
    plt.close();plt.clf()
    plot_full_data(new_fig=True,ref_point=(lats[i],lons[i]), ref_dist=target_distance,time_highlight=date_str, high_no=i,value="temp")
    save_file="Temp_cut_dist_{:02}_{:03}.png".format(target_distance,i)
    plt.savefig(save_file,dpi=120)
    plt.close();plt.clf()
    plot_full_data(new_fig=True,ref_point=(lats[i],lons[i]), ref_dist=target_distance,time_highlight=date_str, high_no=i,value="oxygen", vmin=0,vmax=50,tmin=735600,tmax=735850)
    save_file="Oxygen_cut_dist_{:02}_{:03}.png".format(target_distance,i)
    plt.savefig(save_file,dpi=120)
"""



print "distance figure"

#plt.clf()
fig=plt.figure(figsize=(15,10))   
plt.plot(time,dist_under[:,10],'g')
plt.plot(time,dist_under[:,20],'y')
plt.plot(time,dist_under[:,40],'r')

plt.fill_between(time,dist_under[:,20],dist_under[:,40],facecolor='r')    
plt.fill_between(time,dist_under[:,10],dist_under[:,20],facecolor='y')    
plt.fill_between(time,time*0.0,dist_under[:,10],facecolor='g')    
locs,labels = plt.xticks()
plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%m-%y'))
plt.setp(labels,rotation=45)
plt.ylim((0,points))
plt.yticks(range(0,points,20))
plt.ylabel("Profiles closer than")
plt.xlabel("Time")
plt.legend(['10 km','20 km','40 km'])
#    plt.plot([time[i],time[i]],[plt.ylim()[0],plt.ylim()[1]],color="#ff0000",linewidth=3)
marker_time=mp.dates.date2num(datetime.datetime.strptime("20140821","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#888888",linewidth=3)
marker_time=mp.dates.date2num(datetime.datetime.strptime("20150805","%Y%m%d"))
plt.plot([marker_time,marker_time],[plt.ylim()[0],plt.ylim()[1]],color="#888888",linewidth=3)

#save_file="Distances_{:03}.png".format(i)
#plt.savefig(save_file,dpi=120)

#plt.figure(figsize=(10,15));plot_full_data(value='oxygen',vmin=0,vmax=50,tmin=735600,tmax=735850)
#plt.figure(figsize=(10,15));plot_full_data(value='salt',vmin=11,vmax=13,tmin=735600,tmax=735850)
#Hot SPOTS
