#!/usr/bin/env python
# coding: utf-8
# @author: Yue Sun, UCB

import lakeshore as ls
from OrensteinLab_git.configuration import config_dict
import time
import numpy as np

lakeshore_model = config_dict['Lakeshore Model']

def set_temperature(temperature, lsobj=None, tolerance=0.05, wait_time=0, max_check=750, output=1, on_off=0, rate=0):
    '''
    sets lakeshore setpoint, waits until temperature is within tolerance of setpoint, and waits for soak time before returning.

    args:
        - temperature(float)

    returns: None
    '''
    # Lakeshore initialization
    lsobj, lsobj_passed = get_lsobj(lsobj)

    temp=float(temperature)
    set_ramp(lsobj, output, on_off, rate)
    set_setpoint(temp, lsobj, output)
    time.sleep(0.1)

    # check for stability
    current_temp = []
    for m in range(max_check):
        current_temp.append(read_temperature(lsobj))
        if m >= 10 and abs(np.mean(current_temp[-3:]) - temp) < tolerance:
            time.sleep(wait_time)
            break
        else:
            time.sleep(1)
    if m==max_check-1:
        print(f'Maximum time exceeded. Temperature: {read_temperature(lsobj)}')

    if lsobj_passed==False:
        close_lakeshore(lsobj)

def read_temperature(lsobj=None):
    '''
    reads temperature from lakeshore controller

    args: None

    returns:
        - temp(float):  read temperature
    '''
    lsobj, lsobj_passed = get_lsobj(lsobj)
    if lakeshore_model=='335':
        temp = float(lsobj.query('KRDG?'))
    elif lakeshore_model=='336':
        temp = float(lsobj.query('KRDG?A'))
    if lsobj_passed==False:
        close_lakeshore(lsobj)
    return temp

def set_setpoint(set_temperature, lsobj=None, output=1):
    lsobj, lsobj_passed = get_lsobj(lsobj)
    lsobj.command("SETP "+str(output)+','+str(float(set_temperature)))
    if lsobj_passed==False:
        close_lakeshore(lsobj)

def read_setpoint(lsobj=None):

    lsobj, lsobj_passed = get_lsobj(lsobj)
    if lakeshore_model=='335':
        setpoint = float(lsobj.query('SETP?'))
    elif lakeshore_model=='336':
        setpoint = float(lsobj.query('SETP?A'))
    if lsobj_passed==False:
        close_lakeshore(lsobj)
    return setpoint

def set_ramp(lsobj=None, output=1, on_off=0, rate=0):
    lsobj, lsobj_passed = get_lsobj(lsobj)
    lsobj.command("RAMP "+str(output)+','+str(int(on_off))+','+str(rate))
    if lsobj_passed==False:
        close_lakeshore(lsobj)

def read_ramp(lsobj=None):

    lsobj, lsobj_passed = get_lsobj(lsobj)
    if lakeshore_model=='335':
        output = lsobj.query('RAMP?').split(',')
    elif lakeshore_model=='336':
        output = lsobj.query('RAMP?A').split(',')
    on_off = bool(int(output[0]))
    rate = float(output[1])
    if lsobj_passed==False:
        close_lakeshore(lsobj)
    return [on_off, rate]

def set_range(range, lsobj=None, output=1):
    '''
    sets lakeshore range (0=off, 1=Low, 2=Med, 3=High)

    args:
        - range(int)

    kwargs:
        - output(int): 1 or 2

    returns: None
    '''

    lsobj, lsobj_passed = get_lsobj(lsobj)
    range=int(range)
    if range not in [0,1,2,3]:
        raise ValueError(f'{range} is not a valid range. Please choose from (0=off, 1=Low, 2=Med, 3=High)')
    lsobj.command(f'RANGE {output},{range}')
    if lsobj_passed==False:
        close_lakeshore(lsobj)

def read_range(lsobj=None):
    lsobj, lsobj_passed = get_lsobj(lsobj)
    range = lsobj.query('Range?')
    if lsobj_passed==False:
        close_lakeshore(lsobj)
    return range

def initialize_lakeshore():
    if lakeshore_model == '335':
        return ls.model_335.Model335(57600)
    elif lakeshore_model == '336':
        return ls.model_336.Model336()
    else:
        raise ValueError('Invalid Lakeshore Model.')

def close_lakeshore(lsobj):
    lsobj.disconnect_usb()

def get_lsobj(lsobj):
    '''
    helper function for getting object both within another script or on the command line.
    '''
    lsobj_passed = True
    if lsobj==None:
        lsobj = initialize_lakeshore()
        lsobj_passed=False
    return lsobj, lsobj_passed
