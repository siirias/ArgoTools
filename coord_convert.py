# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 10:27:24 2020

@author: siirias
"""

def coord_string_to_min(coordinates):
    # convert degree and fractions to degree, minutes, fractions.
    # separate with space or comma
    c =  re.sub(","," ", coordinates.strip())
    c = c.split(' ')
    c = [i for i in c if len(i)>0] # skip empties
    result = ""
    for i in c:
        print(i)
        deg = int(float(i))
        minute = 60.0*(float(i) -float(deg))
        result+="{}° {:0.4}' ".format(deg,minute)
    return result