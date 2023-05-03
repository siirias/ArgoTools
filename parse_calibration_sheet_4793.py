# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 18:16:44 2022

Plot some polynomials, for calibration work
@author: siirias

"""
import numpy as np
import matplotlib.pyplot as plt
import re

xmin = 0.0
ymin = 7.0
data =\
"""
SENSOR SERIAL NUMBER: 4793
CALIBRATION DATE: 27-Mar-13
SBE 41cp CONDUCTIVIT Y CALIBRATION DATA
PSS 1978: C(35,15 ,0) = 4.2914 Sieme n s/meter
COEFFICIENTS:
g = -9.821290e - 0 01
h = 1.423607e - 0 01
i = -3.273814e - 0 04
j = 4.338152e - 0 05
CPcor = -9.570 0 e -008
CTcor = 3.250 0 e -006
WBOTC = 6.332 2 e -008
BATH TEMP BATH SAL BATH COND INST FREQ INST COND R ESIDUAL
(ITS-90) (PSU) (Siemens/m) (Hz) (Siemens/m) (Siemens/m)
22.0000 0.0000 0.00000 2631.76 0.00000 0.00000
1.0015 34.8984 2.98244 5286.84 2.98245 0.00001
4.5000 34.8782 3.28999 5487.54 3.28998 -0.00001
15.0000 34.8347 4.27363 6084.44 4.27363 -0.00000
18.5000 34.8250 4.61940 6280.52 4.61940 0.00000
23.9940 34.8142 5.17773 6584.51 5.17772 -0.00000
29.0000 34.8075 5.70104 6856.88 5.70105 0.00001
32.5001 34.8029 6.07393 7044.32 6.07392 -0.00001
"""

data2 =\
"""
SENSOR SERIAL NUMBER: 4793
CALIBRATION DATE: 29-Mar-15
SBE 41cp CONDUCTIVITY CALIBRATION DATA
PSS 1978: C(35,15,0) = 4.2914 Siemens/meter
COEFFICIENTS:
g = -9.821266e-001
h = 1.423310e-001
i = -3.134470e-004
j = 4.250802e-005
CPcor = -9.5700e-008
CTcor = 3.2500e-006
WBOTC = 6.3322e-008
BATH TEMP
(ITS-90)
BATH SAL
(PSU)
BATH COND
(Siemens/m)
INST FREQ
(Hz)
INST COND
(Siemens/m)
RESIDUAL
(Siemens/m)
22.0000
1.0000
4.5000
15.0000
18.5000
23.9940
29.0000
32.5000
0.0000
34.8427
34.8232
34.7813
34.7727
34.7632
34.7579
34.7551
0.00000
2.97801
3.28531
4.26777
4.61321
5.17098
5.69383
6.06652
2631.75
5283.52
5484.15
6080.58
6276.54
6580.37
6852.61
7040.03
0.00000
2.97801
3.28531
4.26778
4.61320
5.17098
5.69384
6.06652
0.00000
0.00000
-0.00000
0.00000
-0.00001
0.00000
0.00001
-0.00001"""



data3 =\
"""SENSOR SERIAL NUMBER: 4793
CALIBRATION DATE: 24-Jan-17
SBE 41cp CONDUCTIVITY CALIBRATION DATA
PSS 1978: C(35,15,0) = 4.2914 Siemens/meter
COEFFICIENTS:
g = -9.836474e-001
h = 1.427120e-001
i = -4.180454e-004
j = 5.025814e-005
CPcor = -9.5700e-008
CTcor = 3.2500e-006
WBOTC = 6.3322e-008
BATH TEMP
(° C)
BATH SAL
(PSU)
BATH COND
(S/m)
INSTRUMENT
OUTPUT (Hz)
INSTRUMENT
COND (S/m)
RESIDUAL
(S/m)
22.0000
1.0000
4.5000
15.0000
18.5000
23.9940
29.0000
32.3568
0.0000
34.7105
34.6908
34.6492
34.6402
34.6305
34.6253
34.6224
0.00000
2.96778
3.27405
4.25328
4.59752
5.15342
5.67455
6.03067
2632.31
5276.89
5477.13
6072.41
6267.98
6571.19
6842.87
7022.28
0.00000
2.96778
3.27405
4.25327
4.59752
5.15342
5.67455
6.03067
0.00000
0.00000
0.00000
-0.00001
-0.00000
0.00001
0.00000
-0.00001
"""
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

d = ParseData(data)
d2 = ParseData(data2)
d3 = ParseData(data3)

max_diff_in_c = 0.0
how_many = 0
print("SENSOR: {}".format(d['ser_num']))
for compare_d in [d,d2,d3]:
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