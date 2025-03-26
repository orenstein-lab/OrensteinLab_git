#!/usr/bin/env python
# coding: utf-8
# @author: Yue Sun, UCB


# from OrensteinLab_git.configuration import CONFIG_DICT
import time
import numpy as np
import sys
import serial
from OrensteinLab_git.configuration import CONFIG_DICT

preamp_port = CONFIG_DICT['Preamp COM Port']

# function to write a command to the instrument
def write_to_instrument(command,obj=None):
    obj, obj_passed = get_obj(obj)
    full_command = command + '\r\n'
    obj.write(full_command.encode()) # Send the command as bytes
    time.sleep(0.1) # small delay to ensure the command is sent
    if obj_passed==False:
        close_preamp(obj)

# function to read a response from the instrument
# def read_from_instrument(obj=None):
#     obj, obj_passed = get_obj(obj)

#     data = obj.readline() # read the response
#     if obj_passed==False:
#         close_preamp(obj)
#     return data.decode().strip() #decode the bytes to string

# def query_instrument(command,obj=None):
#     obj, obj_passed = get_obj(obj)
#     full_command = command + '\r\n'
#     obj.write(full_command.encode()) # Send the command as bytes
#     time.sleep(0.1) # small delay to ensure the command is sent
#     data = obj.readline() # read the response

#     if obj_passed==False:
#         close_preamp(obj)
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
       

def close_preamp(obj):
    obj.close()

def get_obj(obj):
    '''
    helper function for getting object both within another script or on the command line.
    '''
    obj_passed = True
    if obj==None:
        obj = initialize_preamp()
        obj_passed=False
    return obj, obj_passed



