import sys
import telnetlib
import pickle
import time
import numpy as np
from OrensteinLab_git.configuration import CONFIG_DICT

host = CONFIG_DICT['Opticool IP']
port = CONFIG_DICT['Opticool Port']

######################
### Core Functions ###
#####################

def set_opticool_temperature(temperature, optc=None, rate=10, mode=0, wait_time=0, check_stability=True):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

    while True:
        try:
            message=f'TEMP {temperature}, {rate}, {mode}'.encode('ascii')
            optc.write(message+('\r\n').encode('ascii'))
            output=optc.read_some().decode('ascii')

            if check_stability:
                info, optc = get_opticool_temp_info(optc)
                temp, status = info
                while status!='Stable':
                    info, optc = get_opticool_temp_info(optc)
                    temp, status = info
                    time.sleep(1)
                break
        except:
            print('failed to set opticool temperature, reinitializing.')
            optc = initialize_opticool()

    time.sleep(wait_time)

    if obj_passed==False:
        close_opticool(optc)
        return None
    else:
        return optc

def get_opticool_temp_info(optc=None):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False


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

    if obj_passed==False:
        close_opticool(optc)
        return [temp, output_status], None
    else:
        return [temp, output_status], optc

def read_opticool_temperature(optc=None):
    while True:
        try:
            info, optc = get_opticool_temp_info(optc)
            break
        except:
            time.sleep(0.01)
            optc = initialize_opticool()
    return info[0], optc

def set_opticool_field(field, optc=None, rate=110, approach=0, mode=0, wait_time=0, check_stability=True):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

    while True:
        try:
            #print('foo')
            message=f'FIELD {field}, {rate}, {approach}, {mode}'.encode('ascii')
            optc.write(message+('\r\n').encode('ascii'))
            output=optc.read_some().decode('ascii')
            #print('bar')

            if check_stability:
                time.sleep(3)
                info, optc = get_opticool_field_info(optc)
                field, status = info
                while status!='Holding (Driven)':
                    info, optc = get_opticool_field_info(optc)
                    field, status = info
                    time.sleep(1)
                break
            else:
                break
        except:
            print('failed to set opticool field, reinitializing.')
            # print(optc)
            optc = initialize_opticool()
            # print(optc)
            pass

    time.sleep(wait_time)

    if obj_passed==False:
        close_opticool(optc)
        return None
    else:
        return optc

def get_opticool_field_info(optc=None):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False


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

    if obj_passed==False:
        close_opticool(optc)
        return [field, output_status], None
    else:
        return [field, output_status], optc

def read_opticool_field(optc=None):
    while True:
        try:
            info, optc = get_opticool_field_info(optc)
            break
        except:
            time.sleep(0.1)
            optc = initialize_opticool()
    return info[0], optc

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
