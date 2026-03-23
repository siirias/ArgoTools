# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 18:16:44 2022

Plot some polynomials, for calibration work
@author: siirias

IN PROGRESS of CONVERSION TO NEWER VERSIONS
"""
import re
import numpy as np
import PyPDF2 as pp
import datetime as dt

in_dir = "C:\\Data\\EARiseQC\\FloatCalibrationData\\all\\"
in_file1 = "SBE 41cp C3503 22Feb17.pdf"
in_file2 = "41CP-3503.pdf"
in_file3 = "SBE 41cp C3503 12Jun18.pdf"

def readPDF(in_dir, in_file):
    pdf_data = pp.PdfReader(open(in_dir + in_file, 'rb'))
    out_lines = []

    for page in pdf_data.pages:
        rows = {}

        def visitor_text(text, cm, tm, font_dict, font_size):
            if not text.strip():
                return
            x = tm[4]
            y = tm[5]

            # round y a bit so fragments on same visual line end up together
            y_key = round(y, 1)
            rows.setdefault(y_key, []).append((x, text))

        page.extract_text(visitor_text=visitor_text)

        # PDF coordinates usually go bottom->top, so reverse=True
        for y in sorted(rows.keys(), reverse=True):
            pieces = sorted(rows[y], key=lambda v: v[0])

            line = ""
            prev_x = None
            for x, text in pieces:
                # insert space if there is a noticeable gap between fragments
                if prev_x is not None and x - prev_x > 20:
                    line += " "
                line += text
                prev_x = x

            out_lines.append(line)

    return "\n".join(out_lines)
    
def ParseData(data_string):
    """
    Parameters
    ----------
    data_string : TYPE string, copied from calibration datasheet
        DESCRIPTION.

    Returns 
    -------
    None.

    """
    right_page_found = False
    reading_numerics = 'notyet' # 'notyet', 'yes', 'done'
    return_value = {}
    ser_num = 0
    CPcor = 0
    CTcor = 0
    WBOTC = 0
    the_date = 0
    coefficients = []
    numbers = [] # if numbers are each in their own line
    numbers_ln = [] # if each row has six values
    for line in data_string.split('\n'):
        if bool(re.match(".*SERIAL NUMBER",line)):
            ser_num = re.search(":\s*(\d+)",line).groups()[0]
        if bool(re.match(".*CONDUCTIVITY CALIBRATION DATA",line)):
            right_page_found = True
        if bool(re.match('.*CALIBRATION\s+DATE:',line)):
            the_date = re.search(".*CALIBRATION\s+DATE:\s*(.+)\s+SBE.*$", line).groups()[0]
            the_date = dt.datetime.strptime(the_date,'%d-%b-%y')
            
        if(right_page_found): #check these only if we are sure page is right.
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
            if reading_numerics != 'done':
                print("LINE",line)
                vals = re.findall(r"[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?", line)
                print(vals)
                if len(vals) == 6:
                    numbers_ln.append(list(map(float, vals)))
                    reading_numerics = 'yes'
                elif reading_numerics == 'yes':
                    reading_numerics = 'done'
            

    tmp = np.array(numbers_ln).transpose()
    return_value['date'] = the_date
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

data = readPDF(in_dir, in_file1)
data2 = readPDF(in_dir, in_file2)
data3 = readPDF(in_dir, in_file3)

d = ParseData(data)
d2 = ParseData(data2)
d3 = ParseData(data3)

max_diff_in_c = 0.0
how_many = 0
print("SENSOR: {}".format(d['ser_num']))
for compare_d in [d,d2,d3]:
    print("{} vs {}".format(d['date'], compare_d['date']))
    average_diff_in_c = 0.0
    for t,f,c in zip(d['bath_t'],d['inst_freq'],d['bath_c']):
        new_c = freo_to_c(compare_d,f,10.0,t)
        print(new_c, c, c - new_c)
        average_diff_in_c += np.abs(c - new_c)
        how_many += 1
        if np.abs(c - new_c)>max_diff_in_c :
            max_diff_in_c = np.abs(c - new_c)
    average_diff_in_c = average_diff_in_c/float(how_many)
    print ("average error: {}".format(average_diff_in_c))

print("Largest difference in conductivity: {}".format(max_diff_in_c))    