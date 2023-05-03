# -*- coding: utf-8 -*-3.
"""
Created on Wed Sep 13 16:11:18 2017

@author: siirias
"""
import sys
sys.path.insert(0,'D:\\svnfmi_merimallit\\qa\\nemo')
import datetime as dt
import calendar
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import netcdf
import argohelper as ah
from matplotlib import colors as mcolors
import xarray as xr

#runfile('D:/ArgoData/plot_full_data.py', wdir='D:/ArgoData')
dir_argo_data = "C:/Data/ArgoData/ArgosForPlot/GotlandDeep/"
dir_CTD_data = "c:/Data/ArgoData/Siiriaetal2017/"
file_CTD_data = "0316544c.csv"
files_Argo = [\
 'argo-profiles-6902014.nc',
 'argo-profiles-6902019.nc',
 'argo-profiles-6902020.nc']

files_Argo = [\
 'argo-profiles-6902014.nc',
]

draw_images=True
from_beginning=True # If wanting to rerun, just changing the methods of parir matching, put to false to save som time.
variables_plotted={'oxy':'Oxygen','tem':'Temperature','sal':'Salinity'}  # Oxygen, Salinity, Temperature
#km/day range 
day_in_km=3.
max_distance_km=60.
how_many_shown=6
#how_many_analyzed=12 #7  #deprecated
analyze_treshold=35.
argo_style='b.-'
ctd_style='r-'

"""
nc Variable	Variable	Units	Description

metavar1	Cruise		
metavar2	Station		
metavar3	Type		
longitude	Longitude	degrees_east	
latitude	Latitude	degrees_north	
metavar4	Bot. Depth	m	
metavar5	Secchi Depth	m	
date_time	Decimal Gregorian Days of the station	days since 2013-01-01 00:00:00 UTC	Relative Gregorian Days with decimal part
			
var1	PRES	db	
var2	TEMP	deg C	
var3	PSAL	psu	
var4	DOXY	ml/l	
var5	PHOS	umol/l	
var6	TPHS	umol/l	
var7	SLCA	umol/l	
var8	NTRA	umol/l	
var9	NTRI	umol/l	
var10	AMON	umol/l	
var11	NTOT	umol/l	
var12	H2SX	umol/l	
var13	PHPH		
var14	ALKY	meq/l	
var15	CPHL	ug/l	
var16	Year (station date)	
"""
#Gotlands deep
#lon_min=17;lat_min=56;lon_max=22;lat_max=59;
lon_min=16.5;lat_min=55.5;lon_max=22.5;lat_max=59.5;
target_lat=57.32; target_lon=20.05; target_rad=1.8*60000000 #rad in km
def sort_by_best(best_hits):
    for i in range(len(best_hits)):
        for j in range(i,len(best_hits)):
            tmp=best_hits[i]
            if(best_hits[i][1]>best_hits[j][1]):
                best_hits[i]=best_hits[j]
                best_hits[j]=tmp
    return best_hits
        
if from_beginning:
    invalid_val=-10000000000.000
    fmk=pd.read_csv(dir_CTD_data + file_CTD_data)
    press=-1.0*fmk[u'PRES [db]'][:]
    longitude=fmk[u'Longitude [degrees_east]'][:]
    latitude=fmk[ u'Latitude [degrees_north]'][:]
    temperature=fmk[u'TEMP [deg C]'][:]
    salinity=fmk[u'PSAL [psu]'][:]
    oxygen=fmk[u'DOXY [ml/l]'][:]
   #purkka, koska happi eri muodossa:
    for i in range(len(oxygen)):
        if(type(oxygen[i])==str and oxygen[i][0]=='<'):
           oxygen[i]=None 
        else:
           oxygen[i]=float(oxygen[i])

    times_s=fmk[u'yyyy-mm-ddThh:mm'][:]
    times=[]
    for i in range(len(times_s)):
        times.append(dt.datetime.strptime(times_s[i],'%Y-%m-%dT%H:%M'))
    ts=[]
    for i in range(len(times)):
        ts.append(calendar.timegm(times[i].timetuple()))
    ctd_data=ah.split_csv_profiles(press,[longitude,latitude,temperature,salinity,ts,oxygen])
   #tehdän se perus aikajana
    timesx=ctd_data[:,0,5]
    times=[]
    for i in range(len(timesx)):
        times.append(dt.datetime.utcfromtimestamp(timesx[i]))  
    distance_mask=[False]*len(times_s)
    for i in range(len(times_s)):
        if(target_rad>ah.distance((target_lat,target_lon),(latitude[i],longitude[i]))):
            distance_mask[i]=True
    #Ja perus paikkajanatkin
    pressure=ctd_data[:,:,0]
    latitude=ctd_data[:,0,2]
    longitude=ctd_data[:,0,1]
    temperature=ctd_data[:,:,3]
    salinity=ctd_data[:,:,4]
    oxygen=ctd_data[:,:,6]
    
    distance_mask=[False]*len(times)
    for i in range(len(times)):
        if(target_rad>ah.distance((target_lat,target_lon),(latitude[i],longitude[i]))):
            distance_mask[i]=True


    #then plot the argo routes
    plot_legends=False
    plot_routes=True
    
    found_floats=[]
                   
    start=mp.dates.datetime.datetime(1000,5,5)
    end=mp.dates.datetime.datetime(3030,5,5)
    
    
    argo_lats_all=np.array([])
    argo_lons_all=np.array([])
    argo_times_all=np.zeros((0), dtype=np.datetime64)
    argo_depth_all=[]
    argo_tem_all=[]
    argo_sal_all=[]
    argo_oxy_all=[]
    argo_label_all=[]
    argo_name_all=[]
    
    for f in files_Argo:
        argo_lats=np.array([])
        argo_lons=np.array([])
        argo_times=np.array([], dtype=np.datetime64)
        argo_depth=[]
        argo_tem=[]
        argo_sal=[]
        argo_oxy=[]
        argo_label=[]
        argo_name=[]
#            a=qa.PointData(f,1,start,end,"argonc")
        a = xr.open_dataset(dir_argo_data + f )
        print ("File {} has {} profiles".format(f,a['LATITUDE'].shape[0]))
        argo_lats = np.array(a['LATITUDE'])
        argo_lons = np.array(a['LONGITUDE'])
        argo_times = np.array(a['JULD'])
        argo_tem = [x for x in np.array(a['TEMP_ADJUSTED'])]
        argo_sal = [x for x in np.array(a['PSAL_ADJUSTED'])]
        argo_depth = [x for x in np.array(a['PRES_ADJUSTED'])]
        argo_oxy = [x for x in np.array(a['DOXY_ADJUSTED'])]
        argo_name.append(a['PLATFORM_NUMBER'][0])
        
        argo_lats_all=np.concatenate((argo_lats_all,argo_lats))
        argo_lons_all=np.concatenate((argo_lons_all,argo_lons))
        argo_times_all=np.concatenate((argo_times_all,argo_times))
        argo_name_all=np.concatenate((argo_name_all,argo_name))
        argo_label_all+=argo_label
        argo_tem_all+=argo_tem
        argo_sal_all+=argo_sal
        argo_depth_all+=argo_depth
        argo_oxy_all+=argo_oxy
        print("pituus {}".format(len(argo_lats_all)))
        
        argo_lats=argo_lons=argo_times=argo_depth=argo_tem=argo_sal=argo_oxy=argo_label=argo_name=None

    #make the distance mask for argo floats
    argo_distance_mask=[False]*len(argo_lats_all)
    for i in range(len(argo_lats_all)):
        if(target_rad>ah.distance((target_lat,target_lon),(argo_lats_all[i],argo_lons_all[i]))):
            argo_distance_mask[i]=True
        
        
print("CTD measurements\t:",sum(distance_mask))
print("Argo measurements\t:",sum(argo_distance_mask))
for i in argo_name_all:
    if not (i in found_floats):
        found_floats.append(i)


#Lets find out the closest Argo measurement, for each CTD measurement:
best_hits=[]
next_best_hits=[]
#for CTD in range(len(times)):
for Argo in range(len(argo_lats_all[:10])):
    closest_index=0
    closest_ctd=0
    closest_dist=-1.
    closest_time=-1.
    closest_dist_km=-1.
    next_best=[closest_index, closest_dist, closest_time,closest_dist_km, closest_ctd,argo_name_all[closest_index]]
#    for Argo in range(len(argo_lats_all)):
    if not np.isnan(argo_oxy_all[Argo][0]): #If this fails this is non-oxygen profile
        for CTD in range(len(times)):
            themask=np.isfinite(oxygen[CTD])
            #if(True or (not np.isnan(oxygen[CTD][1])) and (not oxygen.mask[CTD][1])):
            if(sum(themask)>1): #Meaning if there is more than one actual oxygen value in the profile.
                Argo_name=argo_name_all[Argo]
                dist_km=ah.distance((argo_lats_all[Argo],argo_lons_all[Argo]),(latitude[CTD],longitude[CTD]))
                dist_days= abs(mp.dates.date2num(argo_times_all[Argo])-mp.dates.date2num(times[CTD]))
                dist=dist_km+day_in_km*dist_days
                if((closest_dist<0. or closest_dist>dist) ): #and dist_km<max_distance_km):
                    next_best=[closest_index, closest_dist, closest_time,closest_dist_km, closest_ctd,argo_name_all[closest_index]]
                    closest_dist=dist
                    closest_dist_km=dist_km
                    closest_time=dist_days
                    closest_index=Argo
                    closest_ctd=CTD
        if(closest_dist<0.):
            print("No hit found for {}! ".format(Argo))
        else:
            best_hits.append([0,0,0,0,0,0])
            best_hits[-1][0]=closest_index
            best_hits[-1][1]=closest_dist
            best_hits[-1][2]=closest_time
            best_hits[-1][3]=closest_dist_km
            best_hits[-1][4]=closest_ctd
            best_hits[-1][5]=argo_name_all[closest_index]
            next_best_hits.append(next_best)
            print("."),
    #else:
        #print "Non oxygen profile: {}".format(Argo)
#        print "for CTD",CTD,"Argo",best_hits[CTD][0],"\t\td(km){:.2f} (d){:.0f}".format(best_hits[CTD][2],best_hits[CTD][3])
print("got from loading!")
        

best_hits=sort_by_best(best_hits)

best_hits_all=best_hits[:]
best_hits=[]
used_CTDs={}
for i in range(len(best_hits_all)):
    if not (best_hits_all[i][4] in used_CTDs):
        used_CTDs[best_hits_all[i][4]]=True
        best_hits.append(best_hits_all[i])


#Let's gather the statistics of wanted pairs:
min_depth_okay=-50.0
how_many_to_skip=0
graph_step=10.0
total_fitness=0.0

#matched_Argos={}
selected_profs={}
pair_statistics=pd.DataFrame(columns=['float',
                                      'number',
                                      'a_lat',
                                      'a_lon',
                                      'c_lat',
                                      'c_lon',
                                      'diff_day',
                                      'diff_km',
                                      'diff_tot',
                                      'a_num',
                                      'c_num'])
tmp=pd.DataFrame([[pd.to_datetime(0),pd.to_datetime(0)]],columns=['a_time','c_time'])
tmp.drop(0) #cludge to empty the system, ensuring the dataformats are right.
pair_statistics=pair_statistics.join(tmp)

for float_no in found_floats:
    found_hits=0
    skipped=0
    name_matches=0
    print("FLOAT:",float_no)
    for No in range(len(best_hits)):
        iargo=best_hits[No][0]
        ictd=best_hits[No][4]
        if(float_no == argo_name_all[iargo] and best_hits[No][1]<analyze_treshold and pressure[ictd].min()<min_depth_okay):
            print("{}: Passed as {} < {}".format(float_no,best_hits[No][1],analyze_treshold))
            name_matches+=1            
            #matched_Argos[iargo]=True
            if((not np.isnan(argo_oxy_all[iargo][0])) and (not np.isnan(oxygen[ictd][0]))):
                if(not skipped>=how_many_to_skip):
                    skipped+=1
                else:
                    skipped=0
                    if float_no not in selected_profs.keys():
                        selected_profs[float_no]={}
                    selected_profs[float_no][found_hits]={'ctd':ictd, 'argo':iargo}
                    N=found_hits
                    found_hits+=1
                    fitness=ah.compare_profiles(pressure[ictd,:],temperature[ictd,:],argo_depth_all[iargo][:],argo_tem_all[iargo][:])
                    total_fitness+=fitness
                    print("\nMatch *{}*".format(N+1))
                    print("profile {} fitness {}, FLOAT:{}".format(No,fitness,argo_name_all[iargo]))
                    print("Distance in km: {} \t Distance in days: {} \t Total: {}".format(best_hits[No][3],best_hits[No][2],best_hits[No][1]))
                    print("Location (CTD) {},{}".format(latitude[ictd],longitude[ictd]))
                    print("Location (Argo) {},{}".format(argo_lats_all[iargo],argo_lons_all[iargo]))
                    print("Time (CTD) {}".format(times[ictd]))
                    print("Time (Argo) {}".format(argo_times_all[iargo]))
                    print("CTD {} ARGO {}\n\n".format(ictd,iargo))
                    pair_statistics=pair_statistics.append({'float':argo_name_all[iargo]},ignore_index=True)
                    idx=len(pair_statistics)-1
                    pair_statistics.loc[idx,'number']=N
                    pair_statistics.loc[idx,'a_lat']=argo_lats_all[iargo]
                    pair_statistics.loc[idx,'a_lon']=argo_lons_all[iargo]
                    pair_statistics.loc[idx,'c_lat']=latitude[ictd]
                    pair_statistics.loc[idx,'c_lon']=longitude[ictd]
                    pair_statistics.loc[idx,'a_time']=pd.to_datetime(str(argo_times_all[iargo]))
                    pair_statistics.loc[idx,'c_time']=pd.to_datetime(str(times[ictd]))
                    pair_statistics.loc[idx,'a_num']=iargo
                    pair_statistics.loc[idx,'c_num']=ictd
                    pair_statistics.loc[idx,'diff_day']=best_hits[No][2]
                    pair_statistics.loc[idx,'diff_km']=best_hits[No][3]
                    pair_statistics.loc[idx,'diff_tot']=best_hits[No][1]
        else:
            if float_no == argo_name_all[iargo]:
                print("{}: Failed as {} !< {}. or {}<{}".format(float_no,best_hits[No][1],analyze_treshold, pressure[ictd].min(),min_depth_okay))


#Sort by argo times
pair_statistics=pair_statistics.sort_values(by='a_time')
#And then plot them:
figure_size=[12,12]
colors=['r','g','b']
for float_no in found_floats:
    fig_float,(ax_t,ax_s,ax_o)=plt.subplots(3,1,figsize=figure_size)
    plt.figure(fig_float.number)
    plt.suptitle("FLOAT {}".format(int(float_no)))
    N=0
    for i,No in pair_statistics.iterrows():
        iargo=No['a_num']
        ictd=No['c_num']
        if(No['float']==float_no and No['number']<how_many_shown): #number to plot the best hits
            ictd=int(No['c_num'])
            iargo=int(No['a_num'])
            for vp in variables_plotted:
                if(variables_plotted[vp]=='Oxygen'):
#                    plt.figure(fig_o.number)
                    themask=np.isfinite(oxygen[ictd])
                    plt.axes(ax_o)
                    plt.plot(oxygen[ictd,themask]+graph_step*N,pressure[ictd,themask],ctd_style,zorder=10)
                    plt.plot(argo_oxy_all[iargo][:]+graph_step*N,-1.0*argo_depth_all[iargo][:],argo_style)
                if(variables_plotted[vp]=='Salinity'):
 #                   plt.figure(fig_s.number)
                    plt.axes(ax_s)
                    plt.plot(salinity[ictd,:]+graph_step*N,pressure[ictd,:],ctd_style,zorder=10)
                    plt.plot(argo_sal_all[iargo][:]+graph_step*N,-1.0*argo_depth_all[iargo][:],argo_style)

                    td_direction=1.0
                    if(No['a_time']>No['c_time']):
                        td_direction=-1.0
                    plt.text(5+graph_step*N,-180,"{:.1f}\n{:.1f} d\n{:.1f} km".format(No['diff_tot'],td_direction*No['diff_day'],No['diff_km']))
                    plt.text(5+graph_step*N,-60,No['a_time'].strftime("%d-%m-%Y"),rotation=90)
                if(variables_plotted[vp]=='Temperature'):
#                    plt.figure(fig_t.number)
                    plt.axes(ax_t)
                    plt.plot(temperature[ictd,:]+graph_step*N,pressure[ictd,:],ctd_style,zorder=10)
                    plt.plot(argo_tem_all[iargo][:]+graph_step*N,-1.0*argo_depth_all[iargo][:],argo_style)
                plt.xlabel(variables_plotted[vp])
                plt.ylabel('Depth')
            N+=1
    plt.savefig("{}_{}.png".format(int(float_no),'matches'))
    plt.savefig("{}_{}.eps".format(int(float_no),'matches'))




#Draw the map locations
# for float_id in found_floats:
#     plt.figure(figsize=[5,5])
#     plt.title("{}".format(float_id))
#     map = Basemap(llcrnrlon=17., llcrnrlat=55., urcrnrlon=23., urcrnrlat=60, resolution='i')
#     map.drawcoastlines()
#     colors=['#0000ff','#00ff00','#ff0000','#000000','#ff00ff','#00ffff','#ffff00','#aaaaff','#aaffaa','#ffaaaa','#aaaaaa','#ffaaff','#aaffff','#ffffaa']
#     colors=mcolors.cnames.keys()
#     lines=[]
#     for i in selected_profs[float_id]:
#         ictd=selected_profs[float_id][i]['ctd']
#         iargo=selected_profs[float_id][i]['argo']
#         ctd_x,ctd_y=map(longitude[ictd],latitude[ictd],[argo_lons_all[iargo],argo_lats_all[iargo]])
#         arg_x,arg_y=map(argo_lons_all[iargo],argo_lats_all[iargo])
#         map.plot([ctd_x,arg_x],[ctd_y,arg_y],marker='',color=colors[i],linestyle='-')
#     plt.legend(range(1,8))

#     for i in selected_profs[float_id]:
#         ictd=selected_profs[float_id][i]['ctd']
#         iargo=selected_profs[float_id][i]['argo']
#         ctd_x,ctd_y=map(longitude[ictd],latitude[ictd],[argo_lons_all[iargo],argo_lats_all[iargo]])
#         arg_x,arg_y=map(argo_lons_all[iargo],argo_lats_all[iargo])
#         map.plot(arg_x,arg_y,marker='.',color=colors[i])
#         map.plot(ctd_x,ctd_y,marker='*',color=colors[i])

#     plt.show()
#     plt.savefig("{}_{}.png".format(int(float_id),"locations"))

#print "With magic number {}(km/day), \taverage fitness: {:.2f} \tavg km diff {:.1f} \tavg d diff {:.2f}".format(day_in_km,total_fitness/how_many_analyzed,avg_km_diff, avg_t_diff)
print("With magic number {}(km/day)".format(day_in_km))



print("Estimating the differences between argo and CTD")

shft=5.0
tot_count=0
colors=['r','g','b']
float_no=0

fig_o,ax_o=plt.subplots(figsize=figure_size)
ax_o.grid(b=True,which='both')
ax_o.set_xticks(np.arange(shft*len(found_floats)))
plt.title("Oxygen differences")

fig_t,ax_t=plt.subplots(figsize=figure_size)
ax_t.grid(b=True,which='both')
ax_t.set_xticks(np.arange(shft*len(found_floats)))
plt.title("Temperature differenses")


fig_s,ax_s=plt.subplots(figsize=figure_size)
ax_s.grid(b=True,which='both')
ax_s.set_xticks(np.arange(shft*len(found_floats)))
plt.title("Salinity differenses")



for float_id in found_floats:
    #average_plots[float_id]=1
    for i in selected_profs[float_id]:
        ictd=selected_profs[float_id][i]['ctd']
        iargo=selected_profs[float_id][i]['argo']
        o,d,orig_o=ah.difference_profile(pressure[ictd],oxygen[ictd],argo_depth_all[iargo],argo_oxy_all[iargo])

        plt.figure(fig_o.number)
        plt.plot(o+shft*tot_count,-1.0*d,'.-',color=colors[float_no])
        
        t,d,orig_t=ah.difference_profile(pressure[ictd],temperature[ictd],argo_depth_all[iargo],argo_tem_all[iargo])
        plt.figure(fig_t.number)
        plt.plot(t+shft*tot_count,-1.0*d,'.-',color=colors[float_no])
        
        s,d,orig_s=ah.difference_profile(pressure[ictd],salinity[ictd],argo_depth_all[iargo],argo_sal_all[iargo])
        plt.figure(fig_s.number)
        plt.plot(s+shft*tot_count,-1.0*d,'.-',color=colors[float_no])
    tot_count+=1
    float_no+=1

plt.figure(fig_o.number)
plt.savefig("{}_diff.png".format(variables_plotted['oxy']))
plt.savefig("{}_diff.eps".format(variables_plotted['oxy']))
plt.figure(fig_t.number)
plt.savefig("{}_diff.png".format(variables_plotted['tem']))
plt.savefig("{}_diff.eps".format(variables_plotted['tem']))
plt.figure(fig_s.number)
plt.savefig("{}_diff.png".format(variables_plotted['sal']))
plt.savefig("{}_diff.eps".format(variables_plotted['sal']))




last_surf=30.0
first_deep=80.0
last_deep=150.0

first_mixed=30.0
last_mixed=80.0


#Plot the scatter Oxygen vs d-oxygen
fig,ax=plt.subplots()
float_no=0
for float_id in found_floats:
    deep_points_sum=0.0
    deep_points=0
    surf_points_sum=0.0
    surf_points=0
    for prof in selected_profs[float_id]:
        ictd=selected_profs[float_id][prof]['ctd']
        iargo=selected_profs[float_id][prof]['argo']
        o,d,orig_o=ah.difference_profile(pressure[ictd],oxygen[ictd],argo_depth_all[iargo],argo_oxy_all[iargo])
        for i in range(len(d)):
            if(i==0 and prof==selected_profs[float_id][0]):
                my_label=float_id
                print("I'm a LABEL!!!")
            else:
                my_label=None
            if(d[i]>last_mixed or d[i]<first_mixed):
                plt.plot(orig_o[i],o[i],'.',color=colors[float_no],label=my_label)
    float_no+=1
    ax.grid(b=True,which='both')
plt.legend()
plt.savefig("scatter_o_vs_do.png")


Analysis=pd.DataFrame(columns=['o_top','o_bot','s_top','s_bot','t_top','t_bot', \
                                'Do_top','Do_bot','Ds_top','Ds_bot','Dt_top','Dt_bot'],index=[])
#count the actual differences
for float_id in found_floats:
    deep_points_sum_o=0.0
    deep_points_sum_o_abs=0.0
    deep_points_o=0
    surf_points_sum_o=0.0
    surf_points_sum_o_abs=0.0
    surf_points_o=0

    deep_points_sum_t=0.0
    deep_points_sum_t_abs=0.0
    deep_points_t=0
    surf_points_sum_t=0.0
    surf_points_sum_t_abs=0.0
    surf_points_t=0

    deep_points_sum_s=0.0
    deep_points_sum_s_abs=0.0
    deep_points_s=0
    surf_points_sum_s=0.0
    surf_points_sum_s_abs=0.0
    surf_points_s=0

    for prof in selected_profs[float_id]:
        ictd=selected_profs[float_id][prof]['ctd']
        iargo=selected_profs[float_id][prof]['argo']
        o,do,orig_o=ah.difference_profile(pressure[ictd],oxygen[ictd],argo_depth_all[iargo],argo_oxy_all[iargo])
        t,dt,orig_t=ah.difference_profile(pressure[ictd],temperature[ictd],argo_depth_all[iargo],argo_tem_all[iargo])
        s,ds,orig_s=ah.difference_profile(pressure[ictd],salinity[ictd],argo_depth_all[iargo],argo_sal_all[iargo])
        for i in range(len(do)):
            if(do[i]<=last_surf):
                surf_points_sum_o+=o[i]
                surf_points_sum_o_abs+=abs(o[i])
                surf_points_o+=1
            if(do[i]>=first_deep and do[i]<=last_deep):
                deep_points_sum_o+=o[i]
                deep_points_sum_o_abs+=abs(o[i])
                deep_points_o+=1
                
        for i in range(len(dt)):
            if(dt[i]<=last_surf):
                surf_points_sum_t+=t[i]
                surf_points_sum_t_abs+=abs(t[i])
                surf_points_t+=1
            if(dt[i]>=first_deep and dt[i]<=last_deep):
                deep_points_sum_t+=t[i]
                deep_points_sum_t_abs+=abs(t[i])
                deep_points_t+=1
                
        for i in range(len(ds)):
            if(ds[i]<=last_surf):
                surf_points_sum_s+=s[i]
                surf_points_sum_s_abs+=abs(s[i])
                surf_points_s+=1
            if(ds[i]>=first_deep and ds[i]<=last_deep):
                deep_points_sum_s+=s[i]
                deep_points_sum_s_abs+=abs(s[i])
                deep_points_s+=1
    Analysis.loc[float_id,['o_top','o_bot']]=[surf_points_sum_o/surf_points_o,deep_points_sum_o/deep_points_o]
    Analysis.loc[float_id,['s_top','s_bot']]=[surf_points_sum_s/surf_points_s,deep_points_sum_s/deep_points_s]
    Analysis.loc[float_id,['t_top','t_bot']]=[surf_points_sum_t/surf_points_t,deep_points_sum_t/deep_points_t]

    Analysis.loc[float_id,['Do_top','Do_bot']]=[surf_points_sum_o_abs/surf_points_o,deep_points_sum_o_abs/deep_points_o]
    Analysis.loc[float_id,['Ds_top','Ds_bot']]=[surf_points_sum_s_abs/surf_points_s,deep_points_sum_s_abs/deep_points_s]
    Analysis.loc[float_id,['Dt_top','Dt_bot']]=[surf_points_sum_t_abs/surf_points_t,deep_points_sum_t_abs/deep_points_t]

    print("\nFor float {}:".format(int(float_id)))
    print("Surface bias (oxyg): {}".format(surf_points_sum_o/surf_points_o))
    print("Surface bias (temperature): {}".format(surf_points_sum_t/surf_points_t)) 
    print("Surface bias (salinity): {}\n".format(surf_points_sum_s/surf_points_s))

    print("bottom bias (oxyg): {}".format(deep_points_sum_o/deep_points_o))
    print("bottom bias (temperature): {}".format(deep_points_sum_t/deep_points_t))
    print("bottom bias (salinity): {}".format(deep_points_sum_s/deep_points_s))
    
    print("Average difference: Surface (oxyg): {}".format(surf_points_sum_o_abs/surf_points_o))
    print("Average difference: Surface (temp): {}".format(surf_points_sum_t_abs/surf_points_t))
    print("Average difference: Surface (sali): {}".format(surf_points_sum_s_abs/surf_points_s))

    print("Average difference: Deep (oxyg): {}".format(deep_points_sum_o_abs/surf_points_o))
    print("Average difference: Deep (temp): {}".format(deep_points_sum_t_abs/surf_points_t))
    print("Average difference: Deep (sali): {}".format(deep_points_sum_s_abs/surf_points_s))

print("\n\nSurface defined as over {} m.".format(last_surf))
print("\n\nbottom defined as under {} m and over {}m.".format(first_deep,last_deep))

print(Analysis[['o_top','s_top','t_top','o_bot','s_bot','t_bot']].to_latex(float_format=lambda x: '{:.2f}'.format(x)))
print(Analysis[['Do_top','Ds_top','Dt_top','Do_bot','Ds_bot','Dt_bot']].to_latex(float_format=lambda x: '{:.2f}'.format(x)))

print("extra analysis")
for i in found_floats:
    print("Float: {}".format(i))
    dat=pair_statistics[pair_statistics['float']==i]
    print("distance max:{:2.3} min:{:.3} avg:{:.3}".format(np.max(dat['diff_km']),np.min(dat['diff_km']),np.average(dat['diff_km'])))
    print("time max:{:.3} min:{:.3} avg:{:.3}".format(np.max(dat['diff_day']),np.min(dat['diff_day']),np.average(dat['diff_day'])))
    print("total max:{:.4} min:{:.3} avg:{:.3}".format(np.max(dat['diff_tot']),np.min(dat['diff_tot']),np.average(dat['diff_tot'])))

print()


plt.figure()
for i in range(len(times)):
    col='b.'
    zorder=0
    themask=np.isfinite(oxygen[i])
    if(sum(themask)>1):
        col='ro'
        zorder=1
        
    plt.plot(longitude[i],latitude[i],col,markerfacecolor='None',zorder=zorder)
plt.plot(argo_lons_all,argo_lats_all,'g',zorder=2)
plt.figure()
for i in range(len(times)):
    col='b.'
    zorder=0
    themask=np.isfinite(oxygen[i])
    if(sum(themask)>1):
        col='ro'
        zorder=1
    plt.plot(times[i],latitude[i],col,markerfacecolor='None',zorder=zorder)

plt.figure()
okprofs=0
for i in range(len(times)):
    col='b.'
    themask=np.isfinite(oxygen[i])
    if(sum(themask)>1):
        okprofs+=1
        plt.plot(oxygen[i][themask]+i,pressure[i][themask],'-')
#few extra numbers:
bh=np.array(best_hits)
for i in found_floats:
    print("Float: {} has {} unique hits".format(i,sum(bh[:,5]==i)))
    num=0
    BY15=0
    Dr=[]
    for j in bh:
        if j[5]==i:
            num+=1
            print("{}\ttot: {:4.2f}\tdist: {:4.2f}\ttime: {:4.2f}\tlat,lon: {:4.2f}-{:4.2f}".format(num,j[1],j[3],j[2],latitude[int(j[4])],longitude[int(j[4])]))
            if ah.distance([target_lat,target_lon],[latitude[int(j[4])],longitude[int(j[4])]])<10.0:
                print("*")
                BY15+=1
            Dr.append(j[1])
            print()
            #print ah.distance([target_lat,target_lon],[latitude[int(j[4])],longitude[int(j[4])]])
    print("\tBY:{}.\t\tDr(6):{:4.2f}\tDr(9):{:4.2f} \tDr(12):{:4.2f}".format(BY15,np.mean(Dr[:6]),np.mean(Dr[:9]),np.mean(Dr[:12])))
    print()
for arg in found_floats:
    print("{} Analyzed profiles: {}".format(arg,len(selected_profs[arg])))

print("FINISHED!")




























