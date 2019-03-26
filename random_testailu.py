# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 18:16:11 2016

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
from time import sleep
"""
distances=plot_full_distance(new_fig=(20,5),save_file="Z_full_distance.png",ref_point=(57.315,20.04))
plt.close()                    

           
plot_full_diff_data(value='temp',col_map='bwr',new_fig=(20,10),save_file="A_temp_diff.png", \
                    distances=distances, vmin=-2,vmax=2)
plt.close()                    

plot_full_data(value='temp',col_map='cool',new_fig=(20,10),save_file="B_temp_full.png")
plt.close()                    


plot_full_diff_data(value='salt',col_map='bwr',new_fig=(20,10),save_file="C_salt_diff.png",\
                    distances=distances, vmin=-1,vmax=1)
plt.close()                    
plot_full_data(value='salt',col_map='cool',new_fig=(20,10),save_file="D_salt_full.png")
plt.close()                    



plot_full_diff_data(value='oxygen',col_map='bwr',new_fig=(20,10),save_file="E_oxygen_diff.png",\
                    distances=distances, vmin=-60,vmax=60)
plt.close()                    
plot_full_data(value='oxygen',col_map='cool',new_fig=(20,10),save_file="F_oxygen_full.png")
plt.close()                    



plot_full_diff_data(value='scatter',col_map='bwr',new_fig=(20,10),save_file="G_scatter_diff.png",\
                    distances=distances, vmin=-0.001,vmax=0.001)
plt.close()                    
plot_full_data(value='scatter',col_map='cool',new_fig=(20,10),save_file="H_scatter_full.png")
plt.close()                    
"""
plot_full_data(value='salt',col_map='cool',new_fig=(20,10),save_file="D_salt_full.png",plot_contour=True)
#plt.close()                    

