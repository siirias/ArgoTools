# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 11:57:10 2020

@author: siirias
"""

import os

files = os.listdir()
for f in files:
    print("{} str_len: {}".format(f,len(f)))
    
