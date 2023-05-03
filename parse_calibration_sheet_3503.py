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
"""SENSOR SERIAL NUMBER: 3503
CALIBRATION DATE: 19-Jun-15
SBE 41cp CONDUCTIVITY CALIBRATION DATA
PSS 1978: C(35,15,0) = 4.2914 Siemens/meter
COEFFICIENTS:
g = -1.028007e+000
h = 1.495265e-001
i = -3.432114e-004
j = 4.724081e-005
CPcor = -9.5700e-008
CTcor = 3.2500e-006
WBOTC = -4.2500e-007
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
4.4999
15.0000
18.5000
23.9940
29.0000
32.5000
0.0000
34.8628
34.8431
34.8003
34.7912
34.7814
34.7763
34.7735
0.00000
2.97956
3.28699
4.26986
4.61540
5.17339
5.69650
6.06937
2627.11
5185.88
5380.60
5959.86
6150.29
6445.61
6710.37
6892.70
0.00000
2.97956
3.28699
4.26986
4.61542
5.17337
5.69648
6.06939
0.00000
-0.00000
-0.00000
0.00001
0.00002
-0.00001
-0.00002
0.00002
"""


data2 =\
"""SENSOR SERIAL NUMBER: 3503
CALIBRATION DATE: 22-Feb-17
SBE 41cp CONDUCTIVITY CALIBRATION DATA
PSS 1978: C(35,15,0) = 4.2914 Siemens/meter
COEFFICIENTS:
g = -1.017400e+000
h = 1.481721e-001
i = -4.263132e-004
j = 5.244764e-005
CPcor = -9.5700e-008
CTcor = 3.2500e-006
WBOTC = -4.2500e-007
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
0.9999
4.5000
15.0001
18.5001
23.9940
29.0000
32.5000
0.0000
34.8153
34.7955
34.7525
34.7432
34.7328
34.7265
34.7218
0.00000
2.97588
3.28295
4.26462
4.60973
5.16696
5.68926
6.06137
2627.11
5205.42
5401.40
5984.32
6175.92
6473.01
6739.24
6922.47
0.00000
2.97589
3.28295
4.26461
4.60973
5.16697
5.68927
6.06136
0.00000
0.00001
-0.00001
-0.00001
0.00000
0.00001
0.00001
-0.00001
"""

data3 =\
"""SENSOR SERIAL NUMBER: 3503
CALIBRATION DATE: 31-Jan-18
SBE 41cp CONDUCTIVITY CALIBRATION DATA
PSS 1978: C(35,15,0) = 4.2914 Siemens/meter
COEFFICIENTS:
g = -1.016771e+000
h = 1.479133e-001
i = -3.481665e-004
j = 4.656400e-005
CPcor = -9.5700e-008
CTcor = 3.2500e-006
WBOTC = -4.2500e-007
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
32.5000
0.0000
34.7454
34.7251
34.6821
34.6727
34.6625
34.6566
34.6531
0.00000
2.97048
3.27697
4.25689
4.60137
5.15765
5.67910
6.05074
2627.14
5201.69
5397.42
5979.64
6171.02
6467.81
6733.84
6916.99
0.00000
2.97048
3.27697
4.25689
4.60137
5.15765
5.67910
6.05074
0.00000
-0.00000
0.00000
-0.00000
0.00000
-0.00000
0.00000
-0.00000
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
        print("{: 10.8f} {: 10.8f} {: 10.8f}".format(new_c, c, c - new_c))
        average_diff_in_c += np.abs(c - new_c)
        how_many += 1
        if np.abs(c - new_c)>max_diff_in_c :
            max_diff_in_c = np.abs(c - new_c)
    average_diff_in_c = average_diff_in_c/float(how_many)
    print ("average error: {}".format(average_diff_in_c))

print("Largest difference in conductivity: {}".format(max_diff_in_c))    