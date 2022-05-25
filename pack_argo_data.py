#! /usr/bin/env python
#Simo Siiria   (simo.siiria at fmi.fi)
#

import numpy
from scipy.io import netcdf as NetCDF
import time
import calendar
import datetime
import pylab
import getopt
import sys
import os
import re
import matplotlib
if(os.environ.get('DISPLAY')==None): #If there is no DISPLAY available, don't try to show figure, just make the file.
    #No display, so set-up matplotlib accordingly
    matplotlib.use('Agg')
import matplotlib.pyplot as plt

#
#Global variables
INVALID_VAL=-10000.0
#column gives the meanings of each row in output file
column={'P':0,'T':1,'S':2,'samples':3,'lon':4,'lat':5,'time':6,'buoy':7,'profile':8,'is_park':9}


#Make_entry should be always used when creating a list which is supposed to go as a line in the output file
def make_entry(P=None,T=None,S=None,lon=None,lat=None,samples=None,time=None,buoy=None,profile=None,is_park=None):
  #make_entry uses column to format a list that fits for the output file.
  result=[0]*10
  result[column['P']]=P
  result[column['T']]=T
  result[column['S']]=S
  result[column['samples']]=samples
  result[column['lon']]=lon
  result[column['lat']]=lat
  result[column['time']]=time
  result[column['buoy']]=buoy
  result[column['profile']]=profile
  result[column['is_park']]=is_park
  for i in range(len(result)):
    if result[i]==None:
      result[i]=INVALID_VAL
  result=map(str,result)
  return result

#
# Main
def main(argv):
  #is_park=0 for points got from profile data, and 
  #is_park=1 for points got from ParkPts data.
  help_text=\
"""
flags:
nothing interesting here yet
"""
  try:                                
    opts, args = getopt.getopt(argv, "h", ["help"])
  except getopt.GetoptError:
    print help_text
    sys.exit(2)                     
  for opts, args in opts:
    if opts in ("-h","--help"):
      print help_text
      sys.exit(2)



  msgfiles=[f for f in os.listdir('.') if re.match('.\d+\.\d+\.msg',f)] 
  f_out=open('profiles.dat','w')
  prof_dat=[]
#
# go through all msg files
  for fname in msgfiles:
    f=open(fname,'r')

    #the file name is asumed to be of form <buoy no>.<profile no.>.msg
    prof_id=re.search('(\d+\.\d+)\.msg',fname).group(1)
    buoy=prof_id.split('.')[0]
    profile=prof_id.split('.')[1]

    msg_data=open(prof_id+'.msg').read()

#Find the time of the profile
#
    try:
      profile_time=re.search('#(.*)Sbe.*NBin',''.join(open(prof_id+'.msg'))).group(1)
      profile_time=calendar.timegm(time.strptime(profile_time.strip(),"%b %d %Y %H:%M:%S"))
      ##Now profile_time is seconds since the epoch GMT.
    except:
      profile_time=INVALID_VAL

#Find the lon lat location of the byou, if available:
#
    try:
      lon=re.search('Fix:\s+([\d\.]+)',''.join(msg_data)).group(1)
      lat=re.search('Fix:\s+[\d\.]+\s+([\d\.]+)',''.join(msg_data)).group(1)
    except:
      lon=str(INVALID_VAL)
      lat=str(INVALID_VAL)

# Read the actual profile data from this file
#
    tmp=read_profile(msg_data)
    for i in tmp:
      prof_dat.append(make_entry( P=i[0],    T=i[1],  S=i[2],        \
                                  samples=i[3], lon=lon,   lat=lat,         \
                                  time=profile_time, buoy=buoy, profile=profile, \
                                  is_park=0))

#Next add the possible Park phase samples:
#
    for i in msg_data.split('\n'):
      try:
        #tmp=[Unix epoc time, mission time, P, T, S]
        #     0             , 1           , 2, 3, 4
        tmp=re.search('ParkPts:.+\s(\d+\s+\d+\s+[\d\.]+\s+[\d\.]+\s+[\d\.]+)[^\d]*\n',i+'\n').group(1).split()
        prof_dat.append(make_entry( P=tmp[2],    T=tmp[3],  S=tmp[4],        \
                                    samples=255.0, lon=lon,   lat=lat,         \
                                    time=tmp[0], buoy=buoy, profile=profile, \
                                    is_park=1))
        #NOTE: Couldn't fin out how many samples are gathered for the park position. Here it is assumed it is
        #      the macimum amount found on profiles (255). This should be ensured tho.
      except:
        pass

  prof_dat=sorted(prof_dat, key=lambda k:k[column['time']])
  prof_dat=sorted(prof_dat, key=lambda k:k[column['profile']])
  prof_dat=sorted(prof_dat, key=lambda k:k[column['buoy']])
  for i in prof_dat:
    f_out.write(" ".join(i)+"\n")
  f.close() 
  f_out.close()
  print "finished!"			



# read_profiles
#
# takes file_data (original msg file read with read()) and returns
# list of lists of profile data [[P T S B] [P T S B] ... ]
# P=pressure, T=temperature, S=salinity, B number of samples in 2 dbar pressure bin forming the P,T,S values

def read_profile(file_data):
#
#first find and read al the lines containing hexadecimal data of the profile
  f=file_data.split('\n')
  first=-1
  for l in range(len(f)):
    if 'NBin' in f[l]:
      hexlines=re.search('NBin\[(\d+)\]',f[l]).groups(1)[0] #grab the number of bins. (with exception of doubles)
      first=l+1
      last=first+int(hexlines)
      #let's check if first is zeros, and maybe multiple zeroes:
      if '00000000000000' in f[first]:
        try:
          multiplier=int(re.search('\[(\d+)\]',f[first]).groups(1)[0])
        except:
          multiplier=1
        first+=1 #at least the first zero line is passed.
        last-=(multiplier-1) #-1 because first occurence is allready removed at start.
  if(first>-1): #something was found
    f=f[first:last]
  else:
    return []
#
#Hexadecimal data read in f, now conversion to human readable form
  data=[]
  for d in f:
    tmp=[0,0,0,0]
    tmp[0]=d[:4]
    tmp[1]=d[4:8]
    tmp[2]=d[8:12]
    tmp[3]=d[12:]
    ##See APEX profiler user manual page 34 (of 52) for reference of conversions
    ##Pressure
    tmp[0]=int(tmp[0],16)
    if(tmp[0]<=0x7FFF):
      tmp[0]=tmp[0]/10.0          #positive
    elif(tmp[0]>=0x8001):
      tmp[0]=(tmp[0]-65536)/10.0  #negative
    elif(tmp[0]==0xFFFE):
      tmp[0]=-0.2                 #near 0
    else:
      tmp[0]=INVALID_VAL          #out of range
        
    ##Temperature
    tmp[1]=int(tmp[1],16)
    if(tmp[1]<=0xEFFF):               #positive
      tmp[1]=tmp[1]/1000.0
    elif(tmp[1]>=0xF001):             #negative
      tmp[1]=(tmp[1]-65536)/1000.0
    elif(tmp[1]==0xFFFE):             #near 0
      tmp[1]=0.002                    
    else:                             #out of range
      tmp[1]=INVALID_VAL              
    ##Salinity
    tmp[2]=int(tmp[2],16)
    if(tmp[2]<=0xEFFF):               #positive
      tmp[2]=tmp[2]/1000.0            
    elif(tmp[2]>=0xF001):             #negative
      tmp[2]=(tmp[2]-65536)/1000.0
    elif(tmp[2]==0xFFFE):             #near 0
      tmp[2]=-0.002
    else:                             #out of range                          
      tmp[2]=INVALID_VAR
    ##Amount of measurements
    tmp[3]=int(tmp[3],16)
    tmp=map(str,tmp)
    data.append(tmp)
  return data



main(sys.argv[1:])
