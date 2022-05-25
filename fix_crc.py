#! /usr/bin/env python
#Simo Siiria   (simo.siiria at fmi.fi)
#
# Usage: fix_crc <filename>
#
#-adds/updates CRC check number after each line of the given file.
#-assumes lines are of format "xxxx(yyy)" and/or "xxxx(yyy) [crcnum]"
#-if old CRC num in [] exists, removes it and updates it.
#-if old CRC num in [] doesn't exist after line, adds it.
#
#-usable for apex floats mission files.
#
#Requires:
#CRCCITT.py by  Cristian NAVALICI cristian.navalici at gmail dot com

import numpy
import time
#import pylab
import getopt
import sys
import os
import re
from CRCCCITT import CRCCCITT


def add_crc(string):
	crc = CRCCCITT('1D0F')  #This is the crc check format used by argo float.
	tmp=string+' ['+hex(crc.calculate(string)).upper()+']\n'
	tmp=re.sub('\[0X','[0x',tmp) #just a small tuning to get the format stay similar to original.
	return tmp
	



def main(argv):
	if (len(argv)!=1):
		print( "usage 'fix_crc.py <filename>'")
		exit(1)
	lines=open(argv[0]).readlines()
	new_lines=[]
	for l in lines:
		tmp=re.sub('\[.*\].*$','',l) #get rid of  possible old crc
		tmp=tmp.rstrip() #and extra spaces and linefeed
		tmp=add_crc(tmp) #then add new crc
		new_lines.append(tmp)
	outp=open(argv[0],'w')
	outp.writelines(new_lines)
main(sys.argv[1:])
