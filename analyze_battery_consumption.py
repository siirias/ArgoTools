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
import scipy.optimize
import cmocean as cmo


input_dir = "D:\\Data\\ArgoData\\ArgoRawData\\"
out_dir = "D:\\Data\\figures\\Battery\\"
show_fitting = False
the_func = lambda x,a,b: (15.0-13.5)*np.exp(b*(x-a))+13.5

show_trigger_voltage = True
presumed_max_voltage = 15.0
trigger_voltage = 13.5

base_colors = [(1.0,0.0,0.0),   (0.0,1.0,0.0),  (0.0,0.0,1.0),
               (0.8,0.2,0.2),   (0.2,0.8,0.2),  (0.2,0.2,0.8),
               (0.99,0.99,0.0),   (0.99,0.0,0.99),  (0.0,0.99,0.99),
               (0.8,0.8,0.2),   (0.8,0.2,0.8),  (0.2,0.8,0.8),
               (0.9,0.5,0.5),   (0.5,0.9,0.5),  (0.5,0.5,0.9),
               (0.4,0.2,0.2),   (0.2,0.4,0.2),  (0.2,0.4,0.2),
               (0.99,0.7,0.7),   (0.7,0.99,0.7),  (0.7,0.7,0.99)
               ]  # 21 Hand picked colors which are more or less easy to tell from eachothers

float_sets = [
        {'n':"f9234_a", 'loc':"f9234\\", 'sensors':'CTD_OB', 'wmo':'6902027','area':'Baltic Proper'},
        {'n':"f9234_b", 'loc':"f9234\\Old_data\\", 'sensors':'CTD_OB', 'wmo':'6902020','area':'Baltic Proper'},
#        {'n':"f9234_c", 'loc':"f9234\\Old_data\\2013_2015\\", 'sensors':'CTD_OB', 'wmo':'6902014','area':'Baltic Proper'},  # never reached the Trigger voltage
        {'n':"f9089_a", 'loc':"f9089\\BAPE2_2014\\", 'sensors':'CTD_OB', 'wmo':'6902018','area':'Bothnian Sea'},
        {'n':"f9089_b", 'loc':"f9089\\BAPE2_2016\\", 'sensors':'CTD_OB', 'wmo':'6902021','area':'Bothnian Sea'},
        {'n':"f9089_c", 'loc':"f9089\\", 'sensors':'CTD_OB', 'wmo':'6902028','area':'Bothnian Sea'},
        {'n':"f8907", 'loc':"f8907\\", 'sensors':'CTD_O', 'wmo':'6903704','area':'N.Baltic Proper'},
        {'n':"f8543", 'loc':"f8543\\", 'sensors':'CTD_O', 'wmo':'6903700','area':'Bay of Bothnia'},
        {'n':"f8540", 'loc':"f8540\\", 'sensors':'CTD_O', 'wmo':'6903701','area':'Baltic Proper'},
#        {'n':"f8348", 'loc':"f8348\\", 'sensors':'CTD_O', 'wmo':'6903695', 'area':'Barents Sea'},  # Barents sea
#        {'n':"f8348_b", 'loc':"f8348\\iceseason18-19\\", 'sensors':'CTD_OT', 'wmo':''},
        {'n':"f7126_a", 'loc':"f7126\\APE1_2012\\", 'sensors':'CTD', 'wmo':'6901901','area':'Bothnian Sea'},
        {'n':"f7126_b", 'loc':"f7126\\APE1_2014-2015\\", 'sensors':'CTD', 'wmo':'6902017','area':'Bothnian Sea'},
        {'n':"f7126_c", 'loc':"f7126\\", 'sensors':'CTD', 'wmo':'6902023','area':'Bothnian Sea'},
#        {'n':"f7087_a", 'loc':"f7087\\APE2_2013\\", 'sensors':'CTD', 'wmo':'6902013','area':'Bothnian Sea'}, # start voltage very low
        {'n':"f7087_b", 'loc':"f7087\\APE2_2016\\", 'sensors':'CTD', 'wmo':'6902022','area':'Bothnian Sea'},
        {'n':"f7087_c", 'loc':"f7087\\APE2_2017\\", 'sensors':'CTD', 'wmo':'6902029','area':'Bothnian Sea'},
        {'n':"f7087_d", 'loc':"f7087\\", 'sensors':'CTD', 'wmo':'6902030','area':'Bothnian Sea'}
        
        ]

print("WMO;Sensors;Software;Start;End;Profiles;Groundings;avg.Depth;Avg.Control a;Area")

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
    bottom_contacts = []
    bc_count = 0
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
                bc_count+=1
            if(voltage and not incomplete):
                voltages.append(voltage)
                cycles.append(int(re.search(".*\.(\d\d\d)\.",file).groups()[0]))
                dates.append(file_time)
                bottom_contacts.append(bc_count)
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
                
        f_s['voltages'] = np.array(voltages)
        f_s['cycles'] = np.array(cycles)
        try:
            f_s['lifetime'] = np.array(list(map(lambda x: (x-dates[0]).days,dates)))
            f_s['dates'] = dates.copy()
        except:
            f_s['lifetime'] = np.array(list(map(lambda x: None,dates)))
            f_s['dates'] = f_s['lifetime'].copy()
        f_s['travelled_depth'] = np.array(travelled_depth)
        f_s['profile_depths'] = np.array(profile_depths)
        f_s['total_control_steps'] = np.array(total_control_steps)
        f_s['control_steps'] = np.array(control_steps)
        f_s['total_control_actions'] = np.array(total_control_actions)
        f_s['control_actions'] = np.array(control_actions)
        f_s['bottom_contacts'] = np.array(bottom_contacts)
    f_s['trigger_index'] = np.argmin(np.abs(f_s['voltages']-trigger_voltage))
    print("{};{};{};{};{};{};{:.1f};{:.1f};{:.1f};{}".format(\
          f_s['wmo'], f_s['sensors'], f_s['type'],\
          dt.datetime.strftime(f_s['dates'][0],'%Y-%m-%d'),\
          dt.datetime.strftime(f_s['dates'][-1],'%Y-%m-%d'),\
          len(f_s['dates']),\
          f_s['bottom_contacts'][-1],\
          np.mean(f_s['profile_depths']),\
          float(f_s['total_control_actions'][-1])/len(f_s['dates']),\
            f_s['area']))

for i in float_sets:
    i['voltages_perc'] = i['voltages']/i['voltages'].max()
    i['fit']= {}
    for xfield in ['travelled_depth', 'lifetime', 'cycles', 'total_control_steps', 'total_control_actions']:
        limit = i['voltages']>13.5
        fitted, covar = scipy.optimize.curve_fit(\
            the_func, \
            i[xfield][limit], i['voltages'][limit],\
            [0.0, -0.001],\
            bounds = ([-np.inf, -0.5],[np.inf, -0.0000001]), maxfev = 10000)
        i['fit'][xfield] = fitted

    
figure_types = [
        {'title':'PercVoltage per depth distance',
         'xlabel':'Travelled depth',
        'ylabel':'Voltage(fraction of maximum)',
        'xfield':'travelled_depth',
        'yfield':'voltages_perc'},

        {'title':'PercVoltage per time',
         'xlabel':'Mission time(days)',
        'ylabel':'Voltage(fraction of maximum)',
        'xfield':'lifetime',
        'yfield':'voltages_perc'},
         
        {'title':'PercVoltage per cycles',
         'xlabel':'Cycle no',
        'ylabel':'Voltage(fraction of maximum)',
        'xfield':'cycles',
        'yfield':'voltages_perc'},

        {'title':'PercVoltage per Control steps',
         'xlabel':'Piston movement steps',
        'ylabel':'Voltage(fraction of maximum)',
        'xfield':'total_control_steps',
        'yfield':'voltages_perc'},

        {'title':'PercVoltage per Control actions',
         'xlabel':'Piston movement actions',
        'ylabel':'Voltage(fraction of maximum)',
        'xfield':'total_control_actions',
        'yfield':'voltages_perc'},
#        ]
#
#figure_types = [
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
        'yfield':'voltages'},

        {'title':'Voltage per Control actions',
         'xlabel':'Piston movement actions',
        'ylabel':'Voltage (V)',
        'xfield':'total_control_actions',
        'yfield':'voltages'}

        ]


#color_scheme = 'by_bottom_contacts'
#color_scheme = 'by_average_profile_depth'
#color_scheme = 'by_average_profile_time'
color_scheme = 'none'

#mark_scheme = 'none'
#mark_scheme = 'area'
mark_scheme = 'sensors'

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
    print(ft['title'])
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
            bc_f = i['bottom_contacts'][-1]/i['cycles'][-1]
            set_color = (bc_f,1.0 - bc_f, 0.5)
        if(color_scheme == 'by_average_profile_depth'):
            bc_f = np.array(i['profile_depths']).mean()/max_profile_depth_average
            set_color = (bc_f, 0.5, 1.0 - bc_f)
        if(color_scheme == 'by_average_profile_time'):
            bc_f = (float(i['lifetime'][-1])/float(i['cycles'][-1]))/max_profile_time_average
            set_color = (0.5, bc_f, 1.0 - bc_f)
        the_mark = None
        if(mark_scheme == 'area'):
            if( i['area'].lower() == 'baltic proper'):
                the_mark = 'x'
            if( i['area'].lower() == 'bothnian sea'):
                the_mark = '.'
            if( i['area'].lower() == 'bay of bothnia'):
                the_mark = 's'
            if( i['area'].lower() == 'n.baltic proper'):
                the_mark = 11
        if(mark_scheme == 'sensors'):
            if( i['sensors'].lower() == 'ctd'):
                the_mark = 'x'
            if( i['sensors'].lower() == 'ctd_o'):
                the_mark = '.'
            if( i['sensors'].lower() == 'ctd_ob'):
                the_mark = 'o'
                
        if('voltages' in i.keys()):
#            the_label = "{} bc:{:.1f} %".format(i['wmo'],100.0*i['bottom_contacts']/i['cycles'][-1])
            the_label = "{}".format(i['wmo'])
            plt.plot(i[ft['xfield']][:],\
                              i[ft['yfield']][:], \
                              label = the_label,\
                              linewidth = lw, color = set_color, \
                              marker = the_mark, markersize = 2.0)
            plt.grid(alpha= 0.25)
            if(show_fitting):
                limit = i['voltages']>13.5
                print(i['fit'][ft['xfield']])
                plt.plot(i[ft['xfield']][limit],\
                                  the_func(np.array(i[ft['xfield']][limit]),\
                                  i['fit'][ft['xfield']][0],\
                                  i['fit'][ft['xfield']][1]), \
                                  linewidth = lw, color = set_color,\
                                  alpha = 0.25, marker = the_mark )
            if(show_trigger_voltage and ft['yfield'] =='voltages_perc'):
                plt.axhline(trigger_voltage/i['voltages'].max(),\
                            color = set_color, alpha = 0.99)
        else:
            print("Different format: {}".format(i['n']))
    if(show_trigger_voltage):
        if(ft['yfield'] =='voltages'):
            plt.axhline(trigger_voltage,color = (0.5,0.5,0.5))
        if(ft['yfield'] =='voltages_perc'):
            plt.axhline(trigger_voltage/presumed_max_voltage,color = (0.5,0.5,0.5))
            
    plt.legend()
    plt.ylabel(ft['ylabel'])
    plt.xlabel(ft['xlabel'])
    plt.title(ft['title'])
    out_filename = "{}.jpg".format(re.sub(" ","_",ft['title']))
    plt.savefig(out_dir+out_filename,\
                        facecolor='w', dpi=300, bbox_inches='tight')

#print table with coefficients, and calulate some derived variables:
# NOTE: variables with averages are calculated based only 
# on part above trigger Vltage
print("WMO;Decay;avg_depth;col/prof;day/prof;cntrolstep/prof")
for i in float_sets:
    i['start_voltage'] = i['voltages'][0]
    i['decay'] = i['fit']['lifetime'][1]
    i['mean_profile_depth'] = i['profile_depths'][:i['trigger_index']].mean()
    i['contact_fraction'] = float(i['bottom_contacts'][i['trigger_index']])/i['cycles'][i['trigger_index']]
    i['mean_days_per_profile'] = i['lifetime'][i['trigger_index']]/float(i['trigger_index'])
    i['mean_control_step_per_profile'] = float(i['total_control_steps'][i['trigger_index']])/i['cycles'][i['trigger_index']]
    i['mean_control_actions_per_profile'] = float(i['total_control_actions'][i['trigger_index']])/i['cycles'][i['trigger_index']]

    print("{};{:.3f};{:.1f};{:.2f};{:.1f};{:.1f}".format(\
          i['wmo'],
          i['decay'], i['mean_profile_depth'],\
          i['contact_fraction'], i['mean_days_per_profile'],\
          i['mean_control_step_per_profile'], \
          i['mean_control_actions_per_profile']))        


scatter_plot_types = [
    {'title':'Mission days (until Vt) depth_days',
     'xlabel':'Average profile depth (m)',
     'ylabel':'Average days between profiles',
     'xfield':'mean_profile_depth',
     'yfield':'mean_days_per_profile',
     'zfield':'lifetime',
     'zlabel':'Days before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.thermal},

    {'title':'Mission days (until Vt) depth_steps',
     'xlabel':'Average profile depth (m)',
     'ylabel':'Average control steps per profile',
     'xfield':'mean_profile_depth',
     'yfield':'mean_control_step_per_profile',
     'zfield':'lifetime',
     'zlabel':'Days before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.thermal},

    {'title':'Mission days (until Vt) depth_contacts',
     'xlabel':'Average profile depth (m)',
     'ylabel':'Fraction of bottom contacts',
     'xfield':'mean_profile_depth',
     'yfield':'contact_fraction',
     'zfield':'lifetime',
     'zlabel':'Days before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.thermal},

    {'title':'Mission days (until Vt) steps_contacts',
     'xlabel':'Average control steps per profile',
     'ylabel':'Fraction of bottom contacts',
     'xfield':'mean_control_step_per_profile',
     'yfield':'contact_fraction',
     'zfield':'lifetime',
     'zlabel':'Days before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.thermal},



    {'title':'Profiles (until Vt) depth_days',
     'xlabel':'Average profile depth (m)',
     'ylabel':'Average days between profiles',
     'xfield':'mean_profile_depth',
     'yfield':'mean_days_per_profile',
     'zfield':'cycles',
     'zlabel':'Profiles before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.solar},

    {'title':'Profiles (until Vt) depth_steps',
     'xlabel':'Average profile depth (m)',
     'ylabel':'Average control steps per profile',
     'xfield':'mean_profile_depth',
     'yfield':'mean_control_step_per_profile',
     'zfield':'cycles',
     'zlabel':'Profiles before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.solar},

    {'title':'Profiles (until Vt) depth_contacts',
     'xlabel':'Average profile depth (m)',
     'ylabel':'Fraction of bottom contacts',
     'xfield':'mean_profile_depth',
     'yfield':'contact_fraction',
     'zfield':'cycles',
     'zlabel':'Profiles before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.solar},

    {'title':'Profiles (until Vt) steps_contacts',
     'xlabel':'Average control steps per profile',
     'ylabel':'Fraction of bottom contacts',
     'xfield':'mean_control_step_per_profile',
     'yfield':'contact_fraction',
     'zfield':'cycles',
     'zlabel':'Profiles before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.solar},



    {'title':'Vertical distance (until Vt) depth_days',
     'xlabel':'Average profile depth (m)',
     'ylabel':'Average days between profiles',
     'xfield':'mean_profile_depth',
     'yfield':'mean_days_per_profile',
     'zfield':'travelled_depth',
     'zlabel':'Meters dived before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.haline},

    {'title':'Vertical distance (until Vt) depth_steps',
     'xlabel':'Average profile depth (m)',
     'ylabel':'Average control steps per profile',
     'xfield':'mean_profile_depth',
     'yfield':'mean_control_step_per_profile',
     'zfield':'travelled_depth',
     'zlabel':'Meters dived before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.haline},

    {'title':'Vertical distance (until Vt) depth_contacts',
     'xlabel':'Average profile depth (m)',
     'ylabel':'Fraction of bottom contacts',
     'xfield':'mean_profile_depth',
     'yfield':'contact_fraction',
     'zfield':'travelled_depth',
     'zlabel':'Meters dived before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.haline},

    {'title':'Vertical distance (until Vt) steps_contacts',
     'xlabel':'Average control steps per profile',
     'ylabel':'Fraction of bottom contacts',
     'xfield':'mean_control_step_per_profile',
     'yfield':'contact_fraction',
     'zfield':'travelled_depth',
     'zlabel':'Meters dived before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.haline},

    {'title':'Vertical distance (until Vt) steps_sV',
     'xlabel':'Average control steps per profile',
     'ylabel':'Starting Voltage(V)',
     'xfield':'mean_control_step_per_profile',
     'yfield':'start_voltage',
     'zfield':'travelled_depth',
     'zlabel':'Meters dived before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.haline},




    {'title':'Control steps (until Vt) profs_sV',
     'xlabel':'Number of profiles',
     'ylabel':'Starting Voltage(V)',
     'xfield':'cycles',
     'yfield':'start_voltage',
     'zfield':'total_control_steps',
     'zlabel':'Control steps before trigger voltage ({})'.format(trigger_voltage),
     'cmap':cmo.cm.haline},

     ]
for ft in scatter_plot_types:
    plt.figure(figsize=(6,5))
    lim_min = float_sets[0][ft['zfield']][i['trigger_index']]
    lim_max = lim_min
    cmap = ft['cmap']
    for i in float_sets:
        if(i[ft['zfield']][i['trigger_index']]>lim_max):
            lim_max = i[ft['zfield']][i['trigger_index']]
        if(i[ft['zfield']][i['trigger_index']]<lim_min):
            lim_min = i[ft['zfield']][i['trigger_index']]
    for i in float_sets:  # mission time as a function of frequency, aveverage depth
        x = i[ft['xfield']]
        y = i[ft['yfield']]
        if(type(x) == np.ndarray):
            x = x[i['trigger_index']]  # if the wanted object is list, not number, get the tirgger value
        if(type(y) == np.ndarray):
            y = y[i['trigger_index']]
        the_color = cmap((i[ft['zfield']][i['trigger_index']]-lim_min)/(lim_max-lim_min))
        plt.plot(x,y,marker = 'o', markersize = 10, color = the_color)
        plt.text(x,y,"{}\n{}\n".format(i['wmo'],i['area']), \
                 horizontalalignment = 'center', fontsize = 6, alpha = 0.5)
    plt.grid(alpha=0.25)
    plt.title(ft['title'])
    plt.xlabel(ft['xlabel'])
    plt.ylabel(ft['ylabel'])
    cb = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=lim_min, vmax=lim_max)))
    cb.set_label(ft['zlabel'])
    out_filename = "{}.jpg".format(re.sub(" ","_",ft['title']))
    plt.savefig(out_dir+out_filename,\
                        facecolor='w', dpi=300, bbox_inches='tight')
