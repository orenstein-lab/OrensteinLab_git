#!/usr/bin/env python
# coding: utf-8
# @author: Yue Sun, UCB

import lakeshore as ls
from OrensteinLab_git.configuration import config_dict, motor_dict, instrument_dict

lakeshore_model = config_dict[]

def read_temperature(lsobj=None):
    '''
    reads temperature from lakeshore controller

    args: None

    returns:
        - temp(float):  read temperature
    '''
    lsobj_passed = True
    if lsobj==None:
        lsobj = initialization_lakeshore()
        lsobj_passed=False
    temp = float(lsobj.query('KRDG?'))
    if lsobj_passed==False:
        ls.close_lakeshore(lsobj)
    return temp

def set_temperature(temperature, lsobj=None, tolerance=0.05, wait_time=0, max_check=750):
    '''
    sets lakeshore setpoint, waits until temperature is within tolerance of setpoint, and waits for soak time before returning.

    args:
        - temperature(float)

    returns: None
    '''
    # Lakeshore initialization
    lsobj_passed=True
    if lsobj==None:
        lsobj = initialize_lakeshore()
        lsobj_passed==False

    temp=float(temperature)
    set_ramp(lsobj, 1, 0, 0)
    set_setpoint(lsobj, 1, temp)
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

def set_lakeshore_range(range, output=1):
    '''
    sets lakeshore range (0=off, 1=Low, 2=Med, 3=High)

    args:
        - range(int)

    kwargs:
        - output(int): 1 or 2

    returns: None
    '''
    range=int(range)
    if range not in [0,1,2,3]:
        raise ValueError(f'{range} is not a valid range. Please choose from (0=off, 1=Low, 2=Med, 3=High)')
    lsobj = ls.initialization_lakeshore335()
    lsobj.command(f'RANGE {output},{range}')
    close_lakeshoress(lsobj)

def read_setpoint(inst):
    return float(inst.query('SETP?'))

def set_setpoint(inst, output, set_temperature):
    inst.command("SETP "+str(output)+','+str(float(set_temperature)))

def read_ramp(inst):
    output = inst.query('RAMP?').split(',')
    on_off = bool(int(output[0]))
    rate = float(output[1])
    return [on_off, rate]

def set_ramp(inst, output, on_off, rate):
    inst.command("RAMP "+str(output)+','+str(int(on_off))+','+str(rate))

def set_range(inst, output, range):
    '''
    sets lakeshore range to low, med, or high

    args:
        - output:   1 or 2
        - range:    0=off, 1=low, 2=med, 3=high
    '''
    inst.command(f'RANGE {output},{range}')

def read_range(inst):
    range = inst.query('Range?')
    return range

def initialize_lakeshore():
    if lakeshore_model == '335':
        return ls.model_335.Model335(57600)
    elif lakeshore_model == '336':
        return ls.model_336.Model336()

def close_lakeshore(lsobj):
    lsobj.disconnect_usb()
