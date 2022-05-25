#! /usr/bin/env python
#Simo Siiria   (simo.siiria at fmi.fi)
#
#Requires:
# check_floats.sh by Simo Siiria  simo.siiria at fmi.fi
#CRCCITT.py by  Cristian NAVALICI cristian.navalici at gmail dot com
# iowtopo2_rev03.nc baltic topografy data. Aquired from IOW (http://www.io-warnemuende.de/topography-of-the-baltic-sea.html)
#

import numpy
from scipy.io import netcdf as NetCDF
import time
import pylab
import getopt
import sys
import os
import re
import subprocess
from subprocess import Popen
from subprocess import PIPE
import matplotlib
if(os.environ.get('DISPLAY')==None): #If there is no DISPLAY available, don't try to show figure, just make the file.
    #No display, so set-up matplotlib accordingly
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from CRCCCITT import CRCCCITT

def sanity_check(parameter,value):
	#returns "" if everything is well based on hand-picked limits for each parameter, if a check is failed, an error string is returned
	#
	# change the check values as appropriate
	if parameter=="ParkPistonPos":
		if int(value)>250:
			return parameter+" > 250"+" ("+value+")"
		if int(value)<1:
			return parameter+" < 1"+" ("+value+")"
	if parameter=="ParkPressure":
		if int(value)>500:
			return parameter+" > 500"+" ("+value+")"
		if int(value)<0:
			return parameter+" < 0"+" ("+value+")"
	if parameter=="DownTime":
		if int(value)>2400:
			return parameter+" > 2400"+" ("+value+")"
		if int(value)<10:
			return parameter+" < 10"+" ("+value+")"
			
	return ""


def add_crc(string):
	crc = CRCCCITT('1D0F')  #This is the crc check format used by argo float.
	tmp=string+' ['+hex(crc.calculate(string)).upper()+']\n'
	tmp=re.sub('\[0X','[0x',tmp) #just a small tuning to get the format stay similar to original.
	return tmp
	



def main(argv):
	input_w=30  #just a parameret telling how long the break-lines will be.
	show_data=False
	modify_cfg=False
	bathymetric_file="iowtopo2_rev03.nc"
	float_name="f7126"
	
	max_over_dive=8.0  #if final ParkPts after dive is this much over the set parking position, the piston position will be adjusted upwards.
	park_points_in_desc=5 #How many first parkPts considered when looking for the maximum depth at descent.(desc_depth)
	max_park_piston_pos=150
	min_park_piston_pos=130
	"""
	There proapbly should be a table/function telling which piston position is good initial value for which desired PsrkingPosition, as that
	determines how deep the float goes initially.  For the moment what I know:
	139 -> 82.0, 82.1, 84.4, 84.0, 83.4, 82.3
	140 -> 80.2, 79.4, 78.3, 78.7, 80.2, 77.5, 77.4, 76.8, 76.2, 78.4, 80.6
	141 -> 77.62, 75.9
	"""
	
	safety_distance=20.0  #in meters, how much difference there must be between set depth, and minimum depth around the float.
	safety_radius=0.06 #in degrees, the radius where minimum/maximum depth is checked (deg change in lat)
	min_meas_depth=0.0	#Minimum depth which can be set as target for the float
	max_meas_depth=90.0	#maximum depth which can be set as target for the float
	base_zoom=0.5 #Initial zoom position
	forced_parameters=[]  #Parameters to be determined manually from command line are given here.
	help_text="""
flags:
-h       , --help        :displays this text
-f <name>, --float        :name of the float. deafult f7126
-s, --show 			:shows the map with the argo
-M, --MODIFY             :changes the mission.cfg file for the argo float. Sets the destination depth based on the bathymetry around the float.
--ParkPressure <num> :forces the PArkPressure to be sot on given number. Requires -M/--MODIFY Flag, can be done also with SetParameter
--PistonPosition <num> :forces the PArkPressure to be sot on given number. Requires -M/--MODIFY Flag, can be done also with SetParameter
--SetParameter <Param=Value> : Sets Param from mission.cfg into Value. Requires -M/--MODIFY Flag.  Example: --SetParameter DownTime=1160
	"""
	try:                                
		opts, args = getopt.getopt(argv, "hf:sM", ["help","float=","show","MODIFY","ParkPressure=","PistonPosition=","SetParameter="])
	except getopt.GetoptError:
		print help_text
		sys.exit(2)                     
	for opts, args in opts:
		if opts in ("-h","--help"):
			print help_text
			sys.exit(2) 
		if opts in ("-f","--float"):
			float_name=args
		if opts in ("-s","--show"):
			if    os.environ.get('DISPLAY')!=None:
				show_data=True
		if opts in ("-M","--MODIFY"):
			modify_cfg=True
		if opts in ("--ParkPressure"):
			forced_parameters.append(("ParkPressure",args.strip()))
		if opts in ("--PistonPosition"):
			forced_parameters.append(("ParkPistonPos",args.strip()))
		if opts in ("--SetParameter"):
			tmp=tuple(args.split('='))
			if(len(tmp)!=2): #Format wasn't X=Y
				print "Illegal SetParameter format. Should be --SetParameter Parameter=Value"
				exit(1)
			forced_parameters.append(tmp)
			
	print '#'*input_w
	print '#'+time.ctime().center(input_w-2)+'#'
	print '#'*input_w
#Aquire the mission.cfg file from seadata.fmi.fi server, and read some data from it
	Popen('scp '+float_name+"@seadata.fmi.fi:mission.cfg .",shell=True)
	mission=open('mission.cfg').readlines()
	park_pressure=float(re.search('ParkPressure\(([\d.]*)',' '.join(mission)).groups(0)[0])

#Lets run the bash script that connects to the seadata.fmi.fi server
	float_data=Popen('./check_floats.sh '+float_name,shell=True,stdout=PIPE,stderr=PIPE).communicate()[0]
	print float_data
	#Next, lets mine from the data the location of the float:
	location=re.search('Fix:[\D]*([\d\.]*)[\D]*([\d\.]*)',float_data).groups()
	location=map(float,location) #change strings to floats
	location_new=[None, None]
	for i in range(len(location)):
		degr=int(location[i])
		min=60.0*(location[i]-float(degr))
		sec=60.0*(min-float(int(min)))
		location_new[i]=[degr,min,sec]
	"""
	loc_string="lon %(lond)d:%(lonm)d:%(lons)0.2f,  lat %(latd)d:%(latm)d:%(lats)0.2f\n" % \
			{"lond":location_new[0][0], "lonm":location_new[0][1], "lons":location_new[0][2], \
			  "latd":location_new[1][0], "latm":location_new[1][1], "lats":location_new[1][2] }
	""" ##This loc_string is deg:min:sec  while we actualy need deg:min (and fractions on min, rather than secs)
	loc_string="lon %(lond)d:%(lonm)0.2f,  lat %(latd)d:%(latm)0.2f\n" % \
			{"lond":location_new[0][0], "lonm":location_new[0][1], \
			  "latd":location_new[1][0], "latm":location_new[1][1] }
	
	print loc_string
#Pick the maximum depth where the float ended
	try:
		park_points=map(float,re.compile('ParkPts:'+'\s+[^\s]+'*6+'\s+([^\s]+)',re.M).findall(float_data))
		max_pp=numpy.max(park_points)
		desc_depth=numpy.max(park_points[0:park_points_in_desc])
	except: #couldn't read the depths for some reason.
		print "WARNING: Unable to read some or all ParkPts"
		park_points=[]
		max_pp=None
		desc_depth=None
#Lets grab the earlier locations, and draw the trajectory:
	#get the data from server:
	trajectory_raw=(Popen('ssh seadata.fmi.fi -l '+float_name+" grep Fix: *.msg",shell=True,stdout=PIPE,stderr=PIPE).communicate()[0]).split('\n')
	trajectory=None
	#format the data on a numpy array:
	for n in trajectory_raw:
		try:
			temp=re.search('Fix:[\D]*([\d\.]*)[\D]*([\d\.]*)',n).groups()
			if trajectory != None :
				trajectory=numpy.vstack( [trajectory, numpy.array([float(temp[0]),float(temp[1])])])
			else:
				trajectory=numpy.array([float(temp[0]),float(temp[1])])
		except:
			pass  #just skip the lines where there are no cordinates
#Let's open the bathymetric file.
	bt_file=NetCDF.NetCDFFile(bathymetric_file,'r')
	bathy=numpy.copy(bt_file.variables['Z_MAX'][:])  #Maximum heights for each point
	x_cord=numpy.copy(bt_file.variables['XT_I'][:])
	y_cord=numpy.copy(bt_file.variables['YT_J'][:])
	bathy[bathy>0.0]=numpy.nan  #Let's remove land points.
#Now let's find the minimum and maximum depth around the float:
	#location=[19.9, 61.5] #DEBUG REMEMBER TO REMOVE
	mg_x,mg_y=numpy.meshgrid(x_cord,y_cord)
	try: #rd is the multiplier for longitude change, so that both lat and lon cover around the same distance.
		rd=1.0/numpy.cos(numpy.radians(location[1]))
	except:
		rd=0.0
		
	proximity_mask=(	(mg_x>location[0]-safety_radius*rd) & (mg_x<location[0]+safety_radius*rd) &\
					(mg_y>location[1]-safety_radius) & (mg_y<location[1]+safety_radius))
	min_depth=bathy[proximity_mask].max()  #max and min vice-versa, because depths are negative
	max_depth=bathy[proximity_mask].min()
	print "Minimum depth around the float: %(Min)0.1f" % {"Min": min_depth}
	print "Maximum depth around the float: %(Max)0.1f" % {"Max": max_depth}
#Pick the proposed measurement depth for the area
	prop_depth=-int((min_depth+safety_distance)/5.0)*5  # This is the proposed set depth, on 5 m steps.
	if(prop_depth>max_meas_depth):
		prop_depth=max_meas_depth
	if(prop_depth<min_meas_depth):
		prop_depth=min_meas_depth

	print "Current set depth: %(cur)0.1f m. Proposed set depth: %(prop)0.1f m" % {"cur":park_pressure, "prop": prop_depth}
	if min_depth > park_pressure-safety_distance :
		print "!"*input_w
		print "FLOAT NEARING TOO SHALLOW AREA!"
		print "Set depth: -%(set)0.1f m, set safety: %(safe)0.1f m, current minimum: %(min)0.1f m !" % {"set":park_pressure,"safe":safety_distance,"min":min_depth}
		print "recommended maximum set depth: %(sd)0.1f m !" % {"sd":prop_depth}
		print "!"*input_w

#Possible tinkering with the mission file comes here:
	for i in range(len(mission)):
		if re.search('ParkPressure',mission[i])!=None:
			tmp='ParkPressure(%1d)' % prop_depth
			mission[i]=add_crc(tmp)
			

		if re.search('ParkPistonPos',mission[i])!=None:
			if park_pressure<desc_depth-max_over_dive: #piston position must be increased:
				print "Float parking ended at %(1)0.2f dbar, while the parking pressure was %(2)d dbar!" % {'1':desc_depth,'2':park_pressure}
				print "suggested to increase piston position!"
				old_pos=float(re.search('\(([\d\.]+)\)',mission[i]).groups()[0])
				new_pos=old_pos+1
				if new_pos>max_park_piston_pos:
					new_pos= max_park_piston_pos
				if new_pos<min_park_piston_pos:
					new_pos= min_park_piston_pos
				tmp='ParkPistonPos(%1d)' % new_pos
				mission[i]=add_crc(tmp)
#Next all the paramteres forced in the command line are checked. Only existing parameteres modified here. No new ones added.
		for param in forced_parameters:
			sanity_msg=sanity_check(param[0],param[1])
			if sanity_msg!="": #Something is not as it should
				print "SANITY VIOLATION!"
				print sanity_msg
				exit(1)
			if re.search(param[0],mission[i])!=None:
				tmp=param[0]+'(%s)' % param[1]
				mission[i]=add_crc(tmp)
				
				
			
	f=open('mission_new.cfg','w')
	f.writelines(mission)
	f.close()

#send the new mission-file for the argo-float
	if(modify_cfg==True):
		print "Setting float ParkPressure from %(1)0.0f dbar into %(2)0.0f dbar." % {'1':park_pressure,'2':prop_depth}
		#backup the earlier mission.cfg
		date_tag=time.strftime("%Y%m%d_%H%M")
		Popen('scp mission.cfg '+float_name+"@seadata.fmi.fi:~/mission_history/mission"+date_tag+".cfg",shell=True)
		#and update it
		Popen('scp mission_new.cfg '+float_name+"@seadata.fmi.fi:mission.cfg",shell=True)
		#check the validity, and siplay the result (wait for a moment to ensure the file is copied properly)
		time.sleep(4)
		Popen('ssh seadata.fmi.fi -l '+float_name+' "chkconfig if=mission.cfg cfg=mission_history/mission'+date_tag+'.cfg"',shell=True)





#Plot the bathymetric data on the figure. Data is cut so that the depths around the safety-debths are emphatized.
	#next two lines cut the bathymetric data on the depts relevant for the float.
	mp=2.0 #multiplier of the safety limit, on plotting.
#	bathy[bathy>-(park_pressure-safety_distance*mp)]=-(park_pressure-safety_distance*mp)
	bathy[bathy<-(park_pressure+safety_distance*mp)]=-(park_pressure+safety_distance*mp)
	if(os.environ.get('DISPLAY')!=None):
		pylab.imshow(bathy,origin='lower',interpolation='nearest',aspect=2.0,extent=[x_cord.min(),x_cord.max(),y_cord.min(),y_cord.max()])
		pylab.colorbar()
		pylab.contour(bathy,origin='lower',interpolation='nearest',aspect=2.0,extent=[x_cord.min(),x_cord.max(),y_cord.min(),y_cord.max()])
		pylab.grid()
		#then plot the trajectory:
		pylab.plot(trajectory[:,0],trajectory[:,1],'g')
		pylab.plot(trajectory[:,0],trajectory[:,1],'w+')
		pylab.plot(location[0],location[1],'ro')
		pylab.plot([location[0] -safety_radius*rd, location[0]+safety_radius*rd],[location[1] -safety_radius, location[1]-safety_radius],'w')
		pylab.plot([location[0] -safety_radius*rd, location[0]+safety_radius*rd],[location[1] +safety_radius, location[1]+safety_radius],'w')
		pylab.plot([location[0] -safety_radius*rd, location[0]-safety_radius*rd],[location[1] -safety_radius, location[1]+safety_radius],'w')
		pylab.plot([location[0]+safety_radius*rd, location[0]+safety_radius*rd],[location[1] -safety_radius, location[1]+safety_radius],'w')
		
		pylab.xlim(location[0]-base_zoom*rd,location[0]+base_zoom*rd )
		pylab.ylim(location[1]-base_zoom,location[1]+base_zoom )
		pylab.title("Argo %(name)s: %(place)s (%(x)0.03f, %(y)0.03f)" % {"name":float_name, "place":loc_string, "x":location[0], "y":location[1]})

		pylab.savefig('latest_loc.png')	
		if show_data:
			pylab.show()
main(sys.argv[1:])