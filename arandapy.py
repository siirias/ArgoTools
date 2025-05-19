# -*- coding: utf-8 -*-
"""
Created on Tue Nov  9 18:06:54 2021

@author: siirias
"""
import re
import datetime as dt
import pandas as pd
import os

def get_aranda_filenames(in_dir, whitelist = None, blacklist = None, cruise_year = '\d\d'):
    in_files_tmp = os.listdir(in_dir)  # all files
    in_files_tmp = [i for i in in_files_tmp if re.match("a{}.*a\.cnv".format(cruise_year), i)] # right types
    # then separate with index
    in_files = []
    for i in in_files_tmp:
        index_no = int(re.search(".*(\d\d\d\d)a\.cnv",i).groups()[0])
        if( whitelist == None or index_no in whitelist): #pick the ones listed
            if(blacklist == None or index_no not in blacklist): #but not in black-listed ones
                in_files.append(i)
    return in_files 
               
def read_aranda_file(file_name):
    lines = open(file_name,'r').readlines()
    end_found = False  #search end to find start of data
    data = []
    columns = []
    long_names = []
    unit_names = []
    station_name = None
    for l in lines:
        #search for the headers
        if(re.match("# name \d?",l)):
            try:
                name = re.search("# name \d?.*=([^:]*)",l).groups()[0].strip()
            except:
                name = None
            try:
                long_name = re.search("# name \d?.*=[^:]*:([^\[]*)",l).groups()[0].strip()
            except:
                long_name = ""
            try:
                unit_name = re.search("# name \d?.*\[(.*)\]",l).groups()[0].strip()
            except:
                unit_name = ""
            columns.append(name)                  
            long_names.append(long_name)                  
            unit_names.append(unit_name)                  
        # search other than column metadata
        if(re.match("\*\* Station name",l)):
            try:
                station_name = re.search("\*\* Station name.*:(.*)",l).groups()[0].strip()
            except:
                print("WARNING: station name failed:",l)
                station_name = "?"

        if(re.match("\*\* Index",l)):
            try:
                station_index = int(re.search("\*\* Index.*:(.*)",l).groups()[0].strip())
            except:
                station_index = 0
            
        
        if(re.match("\*\* Latitude",l)):
            try:
                latitude = re.search("\*\* Latitude.*:(.*)",l).groups()[0].strip()
                latitude = float(re.search("(\d*) \d",latitude).groups()[0]) +\
                           float(re.search("\d* ([\d\.]*)",latitude).groups()[0])/60.0 
            except:
                latitude = 0.0
        if(re.match("\*\* Longitude",l)):
            try:
                longitude = re.search("\*\* Longitude.*:(.*)",l).groups()[0].strip()
                longitude = float(re.search("(\d*) \d",longitude).groups()[0]) +\
                           float(re.search("\d* ([\d\.]*)",longitude).groups()[0])/60.0 
            except:
                longitude = 0.0

        if(re.match("# start_time",l)):  
            try:
                the_time = re.search(\
                         "# start_time = ([a-zA-z]* \d* \d* \d*:\d*:\d*)"\
                         ,l).groups()[0].strip()
                #esim: Oct 15 2020 16:43:44
                the_time = dt.datetime.strptime(the_time,"%b %d %Y %H:%M:%S")
            except:
                print("WARNING, Can't parse time!: {}".format(l))
                the_time = dt.datetime(2000,1,1)

        if(end_found):
            l_t = re.sub("\s\s*"," ",l.strip()).split(" ")
            l_t = list(map(lambda x: float(x),l_t))
            l_t.append(latitude)
            l_t.append(longitude)
            l_t.append(the_time)
            data.append(l_t)
        if(re.match('\*END\*',l)):
            end_found = True
    if(not station_name):
        station_name = "unknown"
        
    columns.append('Lat')
    columns.append('Lon')
    columns.append('Time')
    long_names.append('Latitude')
    long_names.append('Longitude')
    long_names.append('Time')
    ctd_data = pd.DataFrame(data,columns = columns)
    return [ctd_data, station_index, station_name,\
            columns, long_names, unit_names]