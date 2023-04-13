# -*- coding: utf-8 -*-
"""
This script contains functions for analyzing conductivity sensor calibration data.

Functions:
- ParseData: Parses a string containing conductivity sensor calibration data and returns a dictionary of values.
- freo_to_c: Calculates conductivity from frequency using calibration data and other parameters.
- analyze_calibration_data: Analyzes a directory of calibration data files for a single sensor.
- extract_text_from_pdf: Extracts all the text from a given PDF file and returns it as a string.

Usage:
- Parse a single calibration data file with ParseData function.
- Use freo_to_c function to calculate conductivity from frequency and calibration data.
- Use analyze_calibration_data function to analyze a directory of calibration data files for a single sensor.

Author: Simo Siiriä, FMI
"""

import os
import re
import numpy as np
import PyPDF2

def is_float(test_string):
    try:
        float(test_string)
        return True
    except:
        return False

def extract_text_from_pdf(filename):
    """
    Extracts text from a PDF file.

    Args:
        filename (str): The name of the PDF file to extract text from.

    Returns:
        str: The text extracted from the PDF file.

    """    
    if filename.endswith('.pdf'):
        with open(filename, 'rb') as file:
            pdf_reader = PyPDF2.PdfFileReader(file)
            text = ''
            for page in range(pdf_reader.getNumPages()): # Loop through all the pages
                text += pdf_reader.getPage(page).extractText()
        return text
    else:
        return ''


def ParseData(data_string):
    """
    Parses calibration data from a string and returns a dictionary with
    the relevant values.

    Parameters
    ----------
    data_string : string
        The calibration data string copied from the datasheet.

    Returns 
    -------
    dict
        A dictionary with the following keys:
        - 'ser_num': The serial number of the sensor
        - 'coefficients': The calibration coefficients (g, h, i, j)
        - 'CPcor': The CP correction factor
        - 'CTcor': The CT correction factor
        - 'WBOTC': The WBOT correction factor
        - 'bath_t': The bath temperatures used during calibration
        - 'bath_s': The bath salinities used during calibration
        - 'bath_c': The measured conductivity values for each bath
        - 'inst_freq': The instrument frequency used for each bath
        - 'inst_c': The instrument conductivity value for each bath
        - 'resid': The residual values for each bath

    """
    return_value = {}
    ser_num = 0
    CPcor = 0
    CTcor = 0
    WBOTC = 0
    coefficients = []
    numbers = [] # if numbers are each in their own line
    numbers_ln = [] # if each row has six values
    for line in data_string.split('\n'):
        if bool(re.match(".*SERIAL NUMBER",line)):
            ser_num = re.search(":\s*(\d+)",line).groups()[0]
        if bool(re.match(".*[ghij] =",line)):
            tmp_coeff = float(re.search(".*[ghij]\s*=\s*([\-+e\d.]+)",\
                            re.sub('\s','',line)).groups()[0])
            coefficients.append(tmp_coeff)
        if bool(re.match(".*CPcor =",line)):
            CPcor = float(re.search(".*CPcor\s*=\s*([e\-+\d.]+)",\
                            re.sub('\s','',line)).groups()[0])
        if bool(re.match(".*CTcor =",line)):
            CTcor = float(re.search(".*CTcor\s*=\s*([e\-+\d.]+)",\
                            re.sub('\s','',line)).groups()[0])
        if bool(re.match(".*WBOTC =",line)):
            WBOTC = float(re.search(".*WBOTC\s*=\s*([e\-+\d.]+)",\
                            re.sub('\s','',line)).groups()[0])
        if bool(re.match("^[-+\d.e]+$",line)):
            numbers.append(float(line))
        elif bool(re.match("^[-+\d.e\s]+$",line)):
            tmp_line = line.split()
            tmp_line = list(map(float,tmp_line))
            numbers_ln.append(tmp_line)

    if(len(numbers)>len(numbers_ln)):
        tmp = np.array(numbers).reshape((6,-1)) # get Bath T, Bath S, Bath C, Inst Freq, Inst C, Resid
    else:
        tmp = np.array(numbers_ln).transpose()
    
    return_value['ser_num'] = ser_num
    return_value['coefficients'] = coefficients
    return_value['CPcor'] = CPcor 
    return_value['CTcor'] = CTcor 
    return_value['WBOTC'] = WBOTC
    return_value['bath_t'] = tmp[0,:]
    return_value['bath_s'] = tmp[1,:]
    return_value['bath_c'] = tmp[2,:]
    return_value['inst_freq'] = tmp[3,:]
    return_value['inst_c'] = tmp[4,:]
    return_value['resid'] = tmp[5,:]
    return return_value
        
def freo_to_c(data, freq, pressure, temp):
    """
    f = INST FREQ * sqrt(1.0 + WBOTC * t) / 1000.0
    Conductivity = (g + hf^2 + if^3 + jf^4) / (1 + dt + ep ) Siemens/meter
    t = temperature[°C)]; p = pressure[decibars]; d = CTcor; e = CPcor;"
    """
    f = freq*np.sqrt(1.0 + data['WBOTC']*temp)/1000.0
    g = data['coefficients'][0]
    h = data['coefficients'][1]
    i = data['coefficients'][2]
    j = data['coefficients'][3]
    c = (g + h*np.power(f,2) + i*np.power(f,3) + j*np.power(f,4))/(1.0+data['CTcor']*temp + data['CPcor']*pressure)
    return c

def analyze_calibration_data(directory):
    max_diff_in_c = 0.0
    how_many = 0
    the_first = True
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
        file_ok = False
        if filename.endswith('.txt'):
            with open(filepath, 'r') as f:
                contents = f.read()
                file_ok = True
        elif filename.endswith('.pdf'):
            contents = extract_text_from_pdf(filepath)
            file_ok = True
            # Need to find part starting with "SENSOR SERIAL NUMBER"
            # having words:  "CONDUCTIVITY CALIBRATION DATA"
            # and can skip parts after "Date, Slope Correction"
            tmp_cont = contents.split('SENSOR SERIAL NUMBER')
            for t in tmp_cont:
                if "CONDUCTIVITY CALIBRATION DATA" in t:
                    contents = 'SENSOR SERIAL NUMBER' + t
            contents = contents[:contents.rfind('Date, Slope Correction')]
        if file_ok:                
            data = ParseData(contents)
            if the_first:
                print("SENSOR: {}".format(data['ser_num']))
                orig_data = data
                the_first = False
            average_diff_in_c = 0.0
            for t,f,c in zip(orig_data['bath_t'],orig_data['inst_freq'],orig_data['bath_c']):
                new_c = freo_to_c(data,f,10.0,t)
                print("{: 10.8f} {: 10.8f} {: 10.8f}".format(new_c, c, c - new_c))                    
                average_diff_in_c += np.abs(c - new_c)
                how_many += 1
                if np.abs(c - new_c)>max_diff_in_c :
                    max_diff_in_c = np.abs(c - new_c)
            average_diff_in_c = average_diff_in_c/float(how_many)
            print ("average error: {}".format(average_diff_in_c))

    print("Largest difference in conductivity: {}".format(max_diff_in_c))                    


# Example usage
analyze_calibration_data('c:/data/argodata/calib_sheet_test/')
#analyze_calibration_data('c:/data/argodata/calib_sheet_test/pdfs/')
