# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 18:24:48 2018

@author: siirias
"""
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import pandas as pd
import ModelQATools as qa
import argohelper as ah

    
    
    


tot_count=0
colors=['r','g','b']
float_no=0
fig_o,ax_o=plt.subplots()
fig_t,ax_t=plt.subplots()
fig_s,ax_s=plt.subplots()
shft=5.0
for float_id in found_floats:
    average_plots[float_id]=1
    for i in selected_profs[float_id]:
        ictd=selected_profs[float_id][i]['ctd']
        iargo=selected_profs[float_id][i]['argo']
        o,d,orig_o=ah.difference_profile(pressure[ictd],oxygen[ictd],argo_depth_all[iargo],argo_oxy_all[iargo])

        plt.figure(fig_o.number)
        plt.plot(o+shft*tot_count,-1.0*d,'.-',color=colors[float_no])
        ax_o.grid(b=True,which='both')
        ax_o.set_xticks(np.arange(shft*len(found_floats)))
        plt.title("Oxygen differences")
        
        t,d,orig_t=ah.difference_profile(pressure[ictd],temperature[ictd],argo_depth_all[iargo],argo_tem_all[iargo])
        plt.figure(fig_t.number)
        plt.plot(t+shft*tot_count,-1.0*d,'.-',color=colors[float_no])
        ax_t.grid(b=True,which='both')
        ax_t.set_xticks(np.arange(shft*len(found_floats)))
        plt.title("Temperature differenses")
        
        s,d,orig_s=ah.difference_profile(pressure[ictd],salinity[ictd],argo_depth_all[iargo],argo_sal_all[iargo])
        plt.figure(fig_s.number)
        plt.plot(s+shft*tot_count,-1.0*d,'.-',color=colors[float_no])
        ax_s.grid(b=True,which='both')
        ax_s.set_xticks(np.arange(shft*len(found_floats)))
        plt.title("Salinity differenses")

    tot_count+=1
    float_no+=1




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
            if(d[i]>last_mixed or d[i]<first_mixed):
                plt.plot(orig_o[i],o[i],'.',color=colors[float_no])
    float_no+=1
    ax.grid(b=True,which='both')



#count the actual differences
for float_id in found_floats:
    deep_points_sum_o=0.0
    deep_points_o=0
    surf_points_sum_o=0.0
    surf_points_o=0

    deep_points_sum_t=0.0
    surf_points_sum_t=0.0

    deep_points_sum_s=0.0
    surf_points_sum_s=0.0

    for prof in selected_profs[float_id]:
        ictd=selected_profs[float_id][prof]['ctd']
        iargo=selected_profs[float_id][prof]['argo']
        o,d,orig_o=ah.difference_profile(pressure[ictd],oxygen[ictd],argo_depth_all[iargo],argo_oxy_all[iargo])
        for i in range(len(d)):
            if(d[i]<=last_surf):
                surf_points_sum_o+=o[i]
                surf_points_sum_t+=t[i]
                surf_points_sum_s+=s[i]
                surf_points+=1
            if(d[i]>=first_deep and d[i]<=last_deep):
                deep_points_sum_o+=o[i]
                deep_points_sum_t+=t[i]
                deep_points_sum_s+=s[i]
                deep_points+=1

    print "\nFor float {}:".format(float_id)
    print "Surface error (oxyg): {}".format(surf_points_sum_o/surf_points)                
    print "Surface error (temperature): {}".format(surf_points_sum_t/surf_points)                
    print "Surface error (salinity): {}\n".format(surf_points_sum_s/surf_points)                

    print "bottom error (oxyg): {}".format(deep_points_sum_o/deep_points)         
    print "bottom error (temperature): {}".format(deep_points_sum_t/deep_points)         
    print "bottom error (salinity): {}".format(deep_points_sum_s/deep_points)         

print "\n\nSurface defined as over {} m.".format(last_surf)
print "\n\nbottom defined as under {} m and over {}m.".format(first_deep,last_deep)
    
       