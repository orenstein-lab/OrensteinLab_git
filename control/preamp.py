#!/usr/bin/env python
# coding: utf-8
# @author: Yue Sun, UCB


# from OrensteinLab_git.configuration import config_dict
import time
import numpy as np
import sys
import serial
from OrensteinLab_git.configuration import config_dict

preamp_port = config_dict['Preamp COM Port']



# function to write a command to the instrument
def write_to_instrument(command,lsobj=None):
    lsobj, lsobj_passed = get_lsobj(lsobj)
    full_command = command + '\r\n'
    lsobj.write(full_command.encode()) # Send the command as bytes
    time.sleep(0.1) # small delay to ensure the command is sent
    if lsobj_passed==False:
        close_preamp(lsobj)

# function to read a response from the instrument
# def read_from_instrument(lsobj=None):
#     lsobj, lsobj_passed = get_lsobj(lsobj)

#     data = lsobj.readline() # read the response
#     if lsobj_passed==False:
#         close_preamp(lsobj)
#     return data.decode().strip() #decode the bytes to string

# def query_instrument(command,lsobj=None):
#     lsobj, lsobj_passed = get_lsobj(lsobj)
#     full_command = command + '\r\n'
#     lsobj.write(full_command.encode()) # Send the command as bytes
#     time.sleep(0.1) # small delay to ensure the command is sent
#     data = lsobj.readline() # read the response

#     if lsobj_passed==False:
#         close_preamp(lsobj)
#     return data.decode().strip() #decode the bytes to string


def set_v_bias(volt): #input bias voltage in (V)
    
    V_mv = int(volt*1000)
    if (V_mv>5000) or (V_mv<-5000):
        ValueError('Bias voltage out of range')
    else:
        command = 'BSLV'+str(V_mv)
        write_to_instrument(command)
        command = 'BSON1'
        write_to_instrument(command)

def bias_off(): #turn off bias voltage
    command = 'BSON0'
    write_to_instrument(command)


def initialize_preamp():
    found_unit = False
    try:
        ser = serial.Serial(
            port = preamp_port,
            baudrate = 9600,
            partity = serial.PARITY_NONE,
            stopbits = serial.STOPBITS_TWO,
            bytesize = serial.EIGHTBITS,
            timeout = 1
        )
        return ser
    except:
        raise ValueError('Cannot Find Preamp')
       

def close_preamp(lsobj):
    lsobj.close()

def get_lsobj(lsobj):
    '''
    helper function for getting object both within another script or on the command line.
    '''
    lsobj_passed = True
    if lsobj==None:
        lsobj = initialize_preamp()
        lsobj_passed=False
    return lsobj, lsobj_passed



