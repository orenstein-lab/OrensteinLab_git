import sys
import telnetlib
import pickle
import time
import numpy as np
from OrensteinLab_git.configuration import config_dict

host = config_dict['Opticool IP']
port = config_dict['Opticool Port']

######################
### Core Functions ###
#####################

def set_opticool_temperature(temperature, optc=None, rate=10, mode=0, wait_time=0):
    
    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

    while True:
        try:
            message=f'TEMP {temperature}, {rate}, {mode}'.encode('ascii')
            optc.write(message+('\r\n').encode('ascii'))
            output=optc.read_some().decode('ascii')

            while get_opticool_temp_info(optc)[1]!='Stable':
                time.sleep(1)
            break
        except:
            print('failed to set opticool temperature, trying agian.')

    if obj_passed==False:
        close_opticool(optc)

    time.sleep(wait_time)

    #return output

def get_opticool_temp_info(optc=None):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

    while True:
        try:
            message=('TEMP?').encode('ascii')
            optc.write(message+('\r\n').encode('ascii'))
            output=optc.read_some().decode('ascii')
            output_value = output.split(',')[1].strip()
            output_status = output.split(',')[3].strip()
            output_status = output_status.replace("\"","")
            if output_value == 'nan':
                temp = output_value
            else:
                temp = float(output_value)
            break
        except:
            print('failed to get opticool tmep info, trying agian.')

    if obj_passed==False:
        close_opticool(optc)

    return temp, output_status

def read_opticool_temperature(optc=None):
    return get_opticool_temp_info(optc)[0]

def set_opticool_field(field, optc=None, rate=110, approach=0, mode=0, wait_time=0):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

    while True:
        try:
            message=f'FIELD {field}, {rate}, {approach}, {mode}'.encode('ascii')
            optc.write(message+('\r\n').encode('ascii'))
            output=optc.read_some().decode('ascii')

            time.sleep(3)
            while get_opticool_field_info(optc)[1]!='Holding (Driven)':
                time.sleep(1)
            break
        except:
            #print('failed to set opticool field, trying again.')
            pass

    if obj_passed==False:
        close_opticool(optc)

    time.sleep(wait_time)

def get_opticool_field_info(optc=None):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

    while True:
        try:
            message=('FIELD?').encode('ascii')
            optc.write(message+('\r\n').encode('ascii'))
            output=optc.read_some().decode('ascii')
            output_value = output.split(',')[1].strip()
            output_status = output.split(',')[3].strip()
            output_status = output_status.replace("\"","")
            if output_value == 'nan':
                field = output_value
            else:
                field = float(output_value)
            break
        except:
            #print('failed to get opticool field info')
            pass

    if obj_passed==False:
        close_opticool(optc)

    return field, output_status

def read_opticool_field(optc=None):
    return get_opticool_field_info(optc)[0]

def initialize_opticool():
    while True:
        try:
            optc=telnetlib.Telnet(host, port, timeout=15)
            optc.read_until(('Connected to QDInstrument Socket Server.\r\n').encode('ascii'))
            break
        except:
            print('failed to initialize opticool, tryping again')
            time.sleep(0.1)

    return optc

def close_opticool(optc):
    message=b'close\r'
    optc.write(message)
    optc.close()