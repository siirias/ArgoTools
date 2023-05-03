# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 17:36:45 2022

@author: siirias
Extra Plots for deliverable.
"""

import datetime as dt
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
import xarray as xr

output_dir = "C:/Data/ArgoData/Figures/"
fig_dpi = 300

WMOs = [6903709, 6903710]
dirs = ["C:/Data/ArgoData/ArgosForPlot/RBR_BalticProper/",
           "C:/Data/ArgoData/ArgosForPlot/RBR_BothnianSea/"]
for WMO, indir in zip(WMOs, dirs):
    area_name = "."
    if("BalticProper" in indir):
        area_name = "Baltic Proper"
    if("BothnianSea" in indir):
        area_name = "Bothnian Sea"
    
    d_arvor = xr.open_dataset(indir + 
        "/GL_PR_PF_{}.nc".format(
            WMO))
    plt.figure(figsize=(8,4));
    plt.plot(d_arvor['TIME'][:], d_arvor['PRES'][:,0],'.')
    plt.gca().set_ylim(0,10)
    plt.xlabel('Time')
    plt.ylabel('Minimum pressure of the profile (dbar)')
    plt.grid('on')
    plt.title('WMO {}, {}'.format(WMO, area_name))
    print("maximum: {:.3}, minimum {:.3}".format(\
        float(np.max(d_arvor['PRES'][:,0])),\
        float(np.min(d_arvor['PRES'][:,0])))   )
    
    filename = "{}_{}".format("Minimum_pressures",WMO)
    plt.savefig(output_dir+filename+'.png' ,\
            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
    print("saved {}".format(output_dir+filename+'.png'))
    
    
