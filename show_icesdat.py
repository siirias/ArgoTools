# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 11:28:09 2017

@author: siirias
"""

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
import re



icesdat_file=open('ICES_CTD_gotland2011-2017/0817106c.csv','r')
icesdat=icesdat_file.readlines()
icesdat_file.close()
for i in range(len(icesdat)):
    
    icesdat[i]=icesdat[i].split(',')