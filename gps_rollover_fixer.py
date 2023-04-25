# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 17:21:10 2023

@author: siirias
"""

import re
import datetime as dt
import sys
import getopt
import shutil
import os
import sys

# Function to get the regex pattern from the time format string
def get_regexp_string(time_format):
    regex = time_format
    regex = re.sub(r'%Y', r'\\d{4}', regex)
    regex = re.sub(r'%m', r'\\d{1,2}', regex)
    regex = re.sub(r'%d', r'\\d{1,2}', regex)
    regex = re.sub(r'%H', r'\\d{2}', regex)
    regex = re.sub(r'%M', r'\\d{2}', regex)
    regex = re.sub(r'%S', r'\\d{2}', regex)
    regex = re.sub(r'%b', r'[A-Za-z]{3}', regex)
    regex = re.sub(r'%B', r'[A-Za-z]{2,15}', regex)
    regex = re.sub(r'%a', r'[A-Za-z]{3}', regex)
    regex = re.sub(r'%A', r'[A-Za-z]{2,15}', regex)
    regex = re.sub(r'%p', r'(AM|PM)', regex)
    return regex


# Function to find all timestamps in the given string that match any of the given time formats
def find_timestamps(string, time_formats):
    timestamps = []         # List to store datetime objects
    found_strings = []      # List to store matched strings
    matching_patterns = []  # List to store the time format strings that match the matched strings
    for time_format in time_formats:
        pattern = get_regexp_string(time_format)
        # Find all substrings in the string that match the pattern
        matches = re.findall(pattern, string)
        # Convert matched strings to datetime objects and store them along with their corresponding strings and formats
        for match in matches:
            try:
                timestamp = dt.datetime.strptime(match, time_format)
                timestamps.append(timestamp)
                found_strings.append(match)
                matching_patterns.append(time_format)
            except ValueError:
                pass
    # Return a list of tuples containing:
    # the matched strings, datetime objects, and matching time format strings
    return list(zip(found_strings, timestamps, matching_patterns))

# Function to fix timestamps for a single line of text
def fix_stamps_for_line(string, time_formats, gps_time_fix):
    # Fix only timestamps older than this date
    fix_dats_older_than = dt.datetime(2010, 1, 1)
    matches = find_timestamps(string, time_formats)
    new_string = string
    # Loop through all matches and fix the timestamps if they are older than fix_dats_older_than
    for m in matches:
        if(m[1] < fix_dats_older_than):
            new_time = m[1] + gps_time_fix
            new_time_str = new_time.strftime(m[2])
            new_string = re.sub(m[0], new_time_str, new_string)
    # Return the updated line of text
    return new_string


# Function to fix timestamps for the entire text
def fix_stamps_for_text(in_data, time_formats=None, gps_time_fix=None):
    """
    Find all substrings in the given string that match any of the given time formats, 
    and adjust the timestamp according to the given gps_time_fix.
    
    Args:
    string: str, the input text containing timestamps that need to be adjusted
    time_formats: list, optional list of strings specifying the formats of timestamps to be searched for
    gps_time_fix: datetime.timedelta object, optional time adjustment value to be added to timestamps
    
    Returns:
    str, the input string with adjusted timestamps
    """
    # Set default values for time_formats and gps_time_fix if not provided
    if(not gps_time_fix):
        gps_time_fix = dt.timedelta(seconds = 7*24*60*60*1024)
    if(not time_formats):
        time_formats = ['%Y-%m-%d %H:%M:%S', '%b %d %Y %H:%M',
                        '%a %b  %d %H:%M:%S %Y', '%m/%d/%Y %H%M%S']
    # Initialize the new text string
    new_text = ""
    # Loop through each line of the input text
    for string in in_data.split('\n'):
        # Call fix_stamps_for_line() to adjust timestamps in the current line
        new_string = fix_stamps_for_line(string, time_formats, gps_time_fix)
        # Append the modified line to the new text string
        new_text += new_string
        # Add a newline character at the end of the line if there isn't one already
        if(new_text[-1] != '\n'):
            new_text += '\n'
    # Return the modified text string
    return new_text



def main(argv):
    inputfile = ''
    backupfile = ''
    opts, args = getopt.getopt(argv,"hi:b:",["ifile=","backup="])
    for opt, arg in opts:
        if opt == '-h':
           print ('gps_rollover_fixer.py -i <inputfile> -o <outputfile>')
           sys.exit()
        elif opt in ("-i", "--ifile"):
           inputfile = arg
        elif opt in ("-b", "--backup"):
           backupfile = arg
    if not os.path.exists(inputfile):
        print(f'file does not exist: {inputfile}')
        sys.exit(1)
    if backupfile != '' and not os.path.exists(backupfile):
        shutil.copy2(inputfile, backupfile)
    else:
        print('backupfile exists, skipping copy.')
    
    
    in_data = ''.join(open(inputfile,'r').readlines())
    new_data = fix_stamps_for_text(in_data)
    if in_data == new_data:  #nothing changed, so no need to write
        print('No dates needing fixing found.')
    else:
        open(inputfile,'w').write(new_data)
        print(f'rewrote {inputfile}')


if __name__ == "__main__":
   main(sys.argv[1:])



"""
in_dir = 'c:/data/argodata/control_tests/'
in_file = '9568.004.msg'
in_data = ''.join(open(in_dir+in_file,'r').readlines())
print(in_data)
print('-'*80)
new_data = fix_stamps_for_text(in_data)
print(new_data)
"""

