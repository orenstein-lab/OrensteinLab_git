#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
import telnetlib

def connect_opticool():
    #HOST = '127.0.0.1'
    HOST='131.243.163.240'
    PORT = '5000'
    telnetObj=telnetlib.Telnet(HOST,PORT,timeout=15)
    telnetObj.read_until(('Connected to QDInstrument Socket Server.\r\n').encode('ascii'))
    return telnetObj

def read_temperature(telnetObj):
    message=('TEMP?').encode('ascii')
    telnetObj.write(message+('\r\n').encode('ascii'))
    output=telnetObj.read_some().decode('ascii')
    output_value = output.split(',')[1].strip()
    output_status = output.split(',')[3].strip()
    output_status = output_status.replace("\"","")
    if output_value == 'nan':
        return [output_value, output_status]
    else:
        return [float(output_value), output_status]

def set_temperature(telnetObj, set_point, rate, mode):
    message=('TEMP '+str(set_point)+', '+str(rate)+', '+str(mode)).encode('ascii')
    telnetObj.write(message+('\r\n').encode('ascii'))
    output=telnetObj.read_some().decode('ascii')
    return output

def read_field(telnetObj):
    message=('FIELD?').encode('ascii')
    telnetObj.write(message+('\r\n').encode('ascii'))
    output=telnetObj.read_some().decode('ascii')
    output_value = output.split(',')[1].strip()
    output_status = output.split(',')[3].strip()
    output_status = output_status.replace("\"","")
    if output_value == 'nan':
        return [output_value, output_status]
    else:
        return [float(output_value), output_status]

def set_field(telnetObj, set_point, rate, mode):
    message=('FIELD '+str(set_point)+', '+str(rate)+', '+str(mode)+', 1').encode('ascii')
    telnetObj.write(message+('\r\n').encode('ascii'))
    output=telnetObj.read_some().decode('ascii')
    return output

def disconnect_opticool(telnetObj):
    telnetObj.close()

