# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:49:08 2020

@author: siirias
"""

import re
import os
import pandas as pd
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import math


input_dir = "D:\\Data\\ArgoData\\ArgoRawData\\"
out_dir = "D:\\Data\\figures\\Battery\\"
base_colors = [(1.0,0.0,0.0),   (0.0,1.0,0.0),  (0.0,0.0,1.0),
               (0.8,0.2,0.2),   (0.2,0.8,0.2),  (0.2,0.2,0.8),
               (0.99,0.99,0.0),   (0.99,0.0,0.99),  (0.0,0.99,0.99),
               (0.8,0.8,0.2),   (0.8,0.2,0.8),  (0.2,0.8,0.8),
               (0.6,0.5,0.5),   (0.5,0.6,0.5),  (0.5,0.5,0.6),
               (0.4,0.2,0.2),   (0.2,0.4,0.2),  (0.2,0.4,0.2),
               (0.99,0.7,0.7),   (0.7,0.99,0.7),  (0.7,0.7,0.99)
               ]  # 21 Hand picked colors which are more or less easy to tell from eachothers

float_sets = [
        {'n':"f9234_a", 'loc':"f9234\\", 'sensors':'CTD_OB', 'wmo':'6902027'},
        {'n':"f9234_b", 'loc':"f9234\\Old_data\\", 'sensors':'CTD_OB', 'wmo':'6902020'},
        {'n':"f9234_c", 'loc':"f9234\\Old_data\\2013_2015\\", 'sensors':'CTD_OB', 'wmo':'6902014'},
        {'n':"f9089_a", 'loc':"f9089\\BAPE2_2014\\", 'sensors':'CTD_OB', 'wmo':'6902018'},
        {'n':"f9089_b", 'loc':"f9089\\BAPE2_2016\\", 'sensors':'CTD_OB', 'wmo':'6902021'},
        {'n':"f9089_c", 'loc':"f9089\\", 'sensors':'CTD_OB', 'wmo':'6902028'},
        {'n':"f8907", 'loc':"f8907\\", 'sensors':'CTD_O', 'wmo':'6903704'},
        {'n':"f8543", 'loc':"f8543\\", 'sensors':'CTD_O', 'wmo':'6903700'},
        {'n':"f8540", 'loc':"f8540\\", 'sensors':'CTD_O', 'wmo':'6903701'},
#        {'n':"f8348", 'loc':"f8348\\", 'sensors':'CTD_O', 'wmo':''},  # Barents sea
#        {'n':"f8348_b", 'loc':"f8348\\iceseason18-19\\", 'sensors':'CTD_OT', 'wmo':''},
        {'n':"f7126_a", 'loc':"f7126\\APE1_2012\\", 'sensors':'CTD', 'wmo':'6901901'},
        {'n':"f7126_b", 'loc':"f7126\\APE1_2014-2015\\", 'sensors':'CTD', 'wmo':'6902017'},
        {'n':"f7126_c", 'loc':"f7126\\", 'sensors':'CTD', 'wmo':'6902023'},
        {'n':"f7087_a", 'loc':"f7087\\APE2_2013\\", 'sensors':'CTD', 'wmo':'6902013'},
        {'n':"f7087_b", 'loc':"f7087\\APE2_2016\\", 'sensors':'CTD', 'wmo':'6902022'},
        {'n':"f7087_c", 'loc':"f7087\\APE2_2017\\", 'sensors':'CTD', 'wmo':'6902029'},
        {'n':"f7087_d", 'loc':"f7087\\", 'sensors':'CTD', 'wmo':'6902030'}
        
        ]
conversion_to_new_piston_step = 16  # New APX11 software has finer steps,
                                    # as such old steps must be multiplied.
#float_sets = [
#        {'n':"f8348", 'loc':"f8348\\"}
#        
#        ]

def get_science_system_log_name(directory, vital_name):
    file_start = re.search("([^.]+\.\d*\.)",vital_name).groups()[0]
    file_start = re.sub("\.","\.",file_start)
    sci_name = [i for i in os.listdir(directory) \
                if re.match(file_start+".*science_log\.csv$",i)]
    try:
        system_name = [i for i in os.listdir(directory) \
                    if re.match(file_start+".*system_log\.txt$",i)]
    except:
        print("failed",vital_name)
        system_name = ['d']
    return [directory+sci_name[0], directory + system_name[0]]

def bottom_contact(log_lines, float_type = 'apex11'):
    stuck_time = 5.0 # hours from which being stuck is calculated
    stuck_limit = 0.05  # dbar in stuck_time hours 
    windowsize = dt.timedelta(hours=stuck_time/2.0)
    contact_found = False
    park_pressure = []
    park_times = []
    # first gather the parking pressures and times
    if(float_type == 'apex11'):
        for l in log_lines:
            if(re.match('.*CTD_P,.*',l)):
                tmp = re.search('CTD_P,([0-9T]*),([-+.0-9]*)',l).groups()
                time_stamp = dt.datetime.strptime(tmp[0],'%Y%m%dT%H%M%S')
                pressure = float(tmp[1])
                park_pressure.append(pressure)
                park_times.append(time_stamp)
        park_pressure = np.array(park_pressure)
        park_times = np.array(park_times)
    if(float_type == 'apex9'):
        for l in log_lines:
            if(re.match('.*ParkPts?:.*',l)):
                tmp = re.search('ParkPts?:\s+([^\s]{3}\s+\d+\s+\d+\s+\d\d:\d\d:\d\d)(\s+[\d\.-]+){3}',l).groups()
                time_stamp = dt.datetime.strptime(tmp[0],'%b %d %Y %H:%M:%S')
                pressure = float(tmp[1].strip())
                park_pressure.append(pressure)
                park_times.append(time_stamp)
        park_pressure = np.array(park_pressure)
        park_times = np.array(park_times)


    # Then check if these match to a ground contact
    for i in range(len(park_times)):
        t_window = abs(park_times-park_times[i])<windowsize
        max_diff = abs(park_pressure[t_window].min() - \
                   park_pressure[t_window].max())
        if(max_diff)<stuck_limit:
            contact_found = True
    
    return contact_found


print("WMO;Sensors;Software;Start;End;Profiles;Groundings;avg.Depth;Avg.Control actions")
for f_s in float_sets:
    current_input_dir = input_dir+f_s['loc']
    f_s['type'] = "apex9"
    files_to_handle = [i for i in os.listdir(current_input_dir) if re.match(".*\.\d\d\d\.log$",i)]
    # Newer floats have different file setup.
    # if these files are not found, let's search the others
    if(len(files_to_handle)==0):
        f_s['type'] = "apex11"
        files_to_handle = [i for i in os.listdir(current_input_dir) if re.match(".*vitals_log.csv$",i)]
    voltages = []
    cycles = []
    dates = []
    travelled_depth = []
    profile_depths = []
    total_control_steps = []
    control_steps = []
    total_control_actions = []
    control_actions = []
    bottom_contacts = 0
    for file in files_to_handle:
        with open(current_input_dir+file,'r') as f:
            file_time = None
            lines = f.readlines()
            voltage = None
            incomplete = False
            deepest = 0.0
            tmp_control_steps = 0.0
            tmp_control_actions = 0
            #Then let's split on which kind of files we handle:
            if(f_s['type'] == "apex9"):
                # add the other file next to this one
                try:
                    lines = lines +  open(current_input_dir+\
                                     re.sub("\.log",".msg",file),'r').readlines()
                except:
                    #failed to read secondary file, so skip this file
                    incomplete = True
                for l in lines:
                    if(re.match(".*CtdPower.*Volts",l)):
                        voltage = float(re.search(".*CtdPower.* ([\d\.]+)Volts",l).groups()[0])
                        
                    if(re.match("\(.*,\s+\d+\s+sec\)",l)):
                        time_str = re.search("\(([^\)\(]*),\s+\d+\s+sec\)",l).groups()[0]
                        the_time = dt.datetime.strptime(time_str,"%b %d %Y %H:%M:%S")
                        if(not file_time):
                            file_time = the_time
                        if(the_time > file_time):
                            file_time = the_time  # take the latest time in file
                            
                    if(re.match("ParkPts?:.*",l)):
                        this_depth = float(
                                re.search("ParkPts?:(\s+[^\s]+){7}",l)\
                                .groups()[0].strip())
                        if(this_depth>deepest):
                            deepest = this_depth
                    
                    if(re.match(".*PistonMoveAbsWTO.*",l)): # find the control steps
                        start_step, end_step = re.search(\
                                    "PistonMoveAbsWTO\(\)\s+(\d+)->(\d+)"\
                                    ,l).groups()
                        tmp_control_steps += conversion_to_new_piston_step*\
                                    np.abs(int(end_step) -int(start_step))
                        tmp_control_actions += 1
                        
            if(f_s['type'] == "apex11"):
                # add the other file next to this one
                try:
                    lines = lines +  open(\
                            get_science_system_log_name(current_input_dir,file)[0]\
                            ,'r').readlines()
                    lines = lines +  open(\
                            get_science_system_log_name(current_input_dir,file)[1]\
                            ,'r').readlines()
                except:
                    #failed to read secondary file, so skip this file
                    incomplete = True
                    print("failed to add lines")
                for l in lines:
                    if(re.match("VITALS_CORE.*",l)): #last line is latests.
                        voltage = float(re.search(\
                        "VITALS_CORE,[^,]*,[^,]*,[^,]*,([^,]*),"\
                        ,l).groups()[0])
                    if(re.match(".*,\d{8}T\d{6},.*",l)):
                        time_str = re.search(",(\d{8}T\d{6}),",l).groups()[0]
                        the_time = dt.datetime.strptime(time_str,"%Y%m%dT%H%M%S")
                        if(not file_time):
                            file_time = the_time
                        if(the_time > file_time):
                            file_time = the_time  # take the latest time in file
                    if(re.match("CTD_P,.*",l)):
                        this_depth = float(
                                re.search("CTD_P,\d+T\d+,([-\d\.]+)",l)\
                                .groups()[0])
                        if(this_depth>deepest):
                            deepest = this_depth

                    if(re.match(".*Adjusting Buoyancy.*",l)):
                        end_step = int(re.search(\
                                   "Adjusting Buoyancy to[^\d]+(\d+)"\
                                   ,l).groups()[0])
                    if(re.match(".*Buoyancy Start Position.*",l)):
                        start_step = int(re.search(\
                                   "Buoyancy Start Position:\s+(\d+)"\
                                   ,l).groups()[0])
                        tmp_control_steps += np.abs(int(end_step) -int(start_step))
                        tmp_control_actions += 1
                        
            if(bottom_contact(lines, f_s['type'])):
                bottom_contacts+=1
            if(voltage and not incomplete):
                voltages.append(voltage)
                cycles.append(int(re.search(".*\.(\d\d\d)\.",file).groups()[0]))
                dates.append(file_time)
                profile_depths.append(deepest)
                if(len(travelled_depth)==0):
                    travelled_depth.append(deepest)
                else:
                    travelled_depth.append(deepest+travelled_depth[-1])
                if(len(total_control_actions)==0):
                    total_control_actions.append(tmp_control_actions)
                else:
                    total_control_actions.append(tmp_control_actions+total_control_actions[-1])
                control_actions.append(tmp_control_actions)
                if(len(total_control_steps)==0):
                    total_control_steps.append(tmp_control_steps)
                else:
                    total_control_steps.append(tmp_control_steps+total_control_steps[-1])
                control_steps.append(tmp_control_steps)
                
        f_s['voltages'] = voltages.copy()
        f_s['cycles'] = cycles.copy()
        try:
            f_s['lifetime'] = list(map(lambda x: (x-dates[0]).days,dates))
            f_s['dates'] = dates.copy()
        except:
            f_s['lifetime'] = list(map(lambda x: None,dates))
            f_s['dates'] = f_s['lifetime'].copy()
        f_s['travelled_depth'] = travelled_depth.copy()
        f_s['profile_depths'] = profile_depths.copy()
        f_s['total_control_steps'] = total_control_steps.copy()
        f_s['control_steps'] = control_steps.copy()
        f_s['total_control_actions'] = total_control_actions.copy()
        f_s['control_actions'] = control_actions.copy()
        f_s['bottom_contacts'] = bottom_contacts
#    print("From {} to {}".format(f_s['dates'][0], f_s['dates'][-1]))
#    print("distance {} m".format(f_s['travelled_depth'][-1]))
    print("{};{};{};{};{};{};{:.1f};{:.1f};{:.1f}".format(\
          f_s['wmo'], f_s['sensors'], f_s['type'],\
          dt.datetime.strftime(f_s['dates'][0],'%Y-%m-%d'),\
          dt.datetime.strftime(f_s['dates'][-1],'%Y-%m-%d'),\
          len(f_s['dates']),\
          f_s['bottom_contacts'],\
          np.mean(np.array(f_s['profile_depths'])),\
          float(f_s['total_control_actions'][-1])/len(f_s['dates'])))
for i in float_sets:
    i['voltage_perc'] = np.array(i['voltages'])/np.array(i['voltages']).max()
        

    
#figure_types = [
#        {'title':'Voltage per depth distance',
#         'xlabel':'Travelled depth',
#        'ylabel':'Voltage(fraction of maximum)',
#        'xfield':'travelled_depth',
#        'yfield':'voltages_perc'},
#
#        {'title':'Voltage per time',
#         'xlabel':'Mission time(days)',
#        'ylabel':'Voltage(fraction of maximum)',
#        'xfield':'lifetime',
#        'yfield':'voltages_perc'},
#         
#        {'title':'Voltage per cycles',
#         'xlabel':'Cycle no',
#        'ylabel':'Voltage(fraction of maximum)',
#        'xfield':'cycles',
#        'yfield':'voltages_perc'},
#
#        {'title':'Voltage per Control steps',
#         'xlabel':'Piston movement steps',
#        'ylabel':'Voltage(fraction of maximum)',
#        'xfield':'total_control_steps',
#        'yfield':'voltages_perc'}
#
#        ]

figure_types = [
        {'title':'Voltage per depth distance',
         'xlabel':'Travelled depth',
        'ylabel':'Voltage (V)',
        'xfield':'travelled_depth',
        'yfield':'voltages'},

        {'title':'Voltage per time',
         'xlabel':'Mission time(days)',
        'ylabel':'Voltage (V)',
        'xfield':'lifetime',
        'yfield':'voltages'},
         
        {'title':'Voltage per cycles',
         'xlabel':'Cycle no',
        'ylabel':'Voltage (V)',
        'xfield':'cycles',
        'yfield':'voltages'},

        {'title':'Voltage per Control steps',
         'xlabel':'Piston movement steps',
        'ylabel':'Voltage (V)',
        'xfield':'total_control_steps',
        'yfield':'voltages'}

        ]


#color_scheme = 'by_bottom_contacts'
#color_scheme = 'by_average_profile_depth'
color_scheme = 'by_average_profile_time'
color_scheme = 'none'
max_profile_depth_average = 0.0
max_profile_time_average = 0.0
# figure out maximum profile depths etc.
for i in float_sets:
    tmp = np.array(i['profile_depths']).mean()
    if(tmp>max_profile_depth_average):
        max_profile_depth_average = tmp
    tmp = i['lifetime'][-1]/i['cycles'][-1]
    if(tmp>max_profile_time_average):
        max_profile_time_average = tmp
        
for ft in figure_types:
    plt.figure(figsize=(12,8))
    next_color = 0        
    for i in float_sets:
        if(i['type'] == 'apex11'):
            lw = 2.0
        else:
            lw = 1.0
        set_color = base_colors[next_color]
        next_color += 1
        if(color_scheme == 'by_bottom_contacts'):
            bc_f = i['bottom_contacts']/i['cycles'][-1]
            set_color = (bc_f,1.0 - bc_f, 0.5)
        if(color_scheme == 'by_average_profile_depth'):
            bc_f = np.array(i['profile_depths']).mean()/max_profile_depth_average
            set_color = (bc_f, 0.5, 1.0 - bc_f)
        if(color_scheme == 'by_average_profile_time'):
            bc_f = i['lifetime'][-1]/i['cycles'][-1]/max_profile_time_average
            set_color = (0.5, bc_f, 1.0 - bc_f)
            
        if('voltages' in i.keys()):
#            the_label = "{} bc:{:.1f} %".format(i['wmo'],100.0*i['bottom_contacts']/i['cycles'][-1])
            the_label = "{}".format(i['wmo'])
            plt.plot(np.array(i[ft['xfield']][:]),\
                              i[ft['yfield']][:], \
                              label = the_label,\
                              linewidth = lw, color = set_color)
        else:
            print("Different format: {}".format(i['n']))
        plt.legend()
        plt.ylabel(ft['ylabel'])
        plt.xlabel(ft['xlabel'])
        plt.title(ft['title'])
        out_filename = "{}.jpg".format(re.sub(" ","_",ft['title']))
        plt.savefig(out_dir+out_filename,\
                            facecolor='w', dpi=300, bbox_inches='tight')

