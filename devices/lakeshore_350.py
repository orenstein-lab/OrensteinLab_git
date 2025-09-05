#!/usr/bin/env python
# coding: utf-8
# @author: Yue Sun, UCB

import lakeshore as ls
from OrensteinLab_git.configuration import CONFIG_DICT
import time
import numpy as np

INPUT = CONFIG_DICT['Lakeshore Input']
OUTPUT = CONFIG_DICT['Lakeshore Output']
INPUTDICT = {'None':0, 'A':1, 'B':2, 'C':3, 'D':4}

def set_temperature(temperature, lsobj=None, tolerance=0.1, avg_time=3, wait_time=0, max_check=750, output=OUTPUT, input=INPUT, on_off=0, rate=0, check_stability=True):
    '''
    sets lakeshore setpoint, waits until temperature is within tolerance of setpoint, and waits for soak time before returning.

    args:
        - temperature(float)

    returns: None
    '''
    # Lakeshore initialization
    lsobj, lsobj_passed = get_lsobj(lsobj)

    temp=float(temperature)
    #lsobj = set_ramp(lsobj, output, on_off, rate)b
    set_setpoint(temp, lsobj, output)
    time.sleep(0.1)

    # check for stability
    if check_stability==True:
        current_temp = []
        for m in range(max_check):
            current_temp.append(read_temperature(lsobj, input)[0])
            if m >= 5*avg_time and abs(np.mean(current_temp[-avg_time:]) - temp) < tolerance:
                time.sleep(wait_time)
                break
            else:
                time.sleep(1)
        if m==max_check-1:
            time.sleep(wait_time)
            print(f'Maximum time exceeded. Temperature: {read_temperature(lsobj, input)[0]}')

    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return None
    else:
        return lsobj

def read_temperature(lsobj=None, input=INPUT):
    '''
    reads temperature from lakeshore controller

    args: None

    returns:
        - temp(float):  read temperature
    '''
    lsobj, lsobj_passed = get_lsobj(lsobj)
    temp = float(lsobj.query(f'KRDG? {input}'))

    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return temp, None
    else:
        return temp, lsobj

def check_lakeshore_stability(lsobj=None, tolerance=0.01, max_check=750, avg_time=30, wait_time=30, input=INPUT):

    current_temp = []
    for m in range(max_check):
        current_temp.append(read_temperature(lsobj, input)[0])
        if m >= 3*avg_time and np.std(current_temp[-avg_time:]) < tolerance:
            break
        else:
            time.sleep(1)

    time.sleep(wait_time)

def set_setpoint(set_temperature, lsobj=None, output=OUTPUT):
    lsobj, lsobj_passed = get_lsobj(lsobj)
    lsobj.command("SETP "+str(output)+','+str(float(set_temperature)))
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return None
    else:
        return lsobj

def read_setpoint(lsobj=None, output=OUTPUT):

    lsobj, lsobj_passed = get_lsobj(lsobj)

    setpoint = float(lsobj.query(f'SETP? {output}'))

    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return setpoint, lsobj
    return setpoint, lsobj

def set_ramp(on_off, rate, lsobj=None, output=OUTPUT):
    '''
    args:
        - output: specifies which output control loop to configure. 1 or 2
        - on_off: specifies whether ramping is 0=off, 1=on
        - rate:   specifies setpoint ramp rate in K/min
    '''
    lsobj, lsobj_passed = get_lsobj(lsobj)
    lsobj.command("RAMP "+str(output)+','+str(int(on_off))+','+str(rate))
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return None
    else:
        return lsobj

def read_ramp(lsobj=None, output=OUTPUT):

    lsobj, lsobj_passed = get_lsobj(lsobj)

    output = lsobj.query(f'RAMP? {output}').split(',')

    on_off = bool(int(output[0]))
    rate = float(output[1])
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return [on_off, rate], None
    return [on_off, rate], lsobj

def set_range(range, lsobj=None, output=OUTPUT):
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
        return None
    else:
        return lsobj

def read_range(lsobj=None, output=OUTPUT):
    lsobj, lsobj_passed = get_lsobj(lsobj)
    range = lsobj.query(f'RANGE? {output}')
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return range, None
    return range, lsobj

def set_pid(p, i, d, lsobj=None, output=OUTPUT):

    lsobj, lsobj_passed = get_lsobj(lsobj)
    lsobj.command(f'PID {output}, {p}, {i}, {d}')
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return None
    else:
        return lsobj

def read_pid(lsobj=None, output=OUTPUT):

    lsobj, lsobj_passed = get_lsobj(lsobj)
    pid = lsobj.query(f'PID? {output}')
    p, i, d = [float(i) for i in pid.split(',')]
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return [p, i, d], None
    return [p, i, d], lsobj

def set_outmode(mode, lsobj=None, output=OUTPUT, input=INPUT, powerup_enable=1):
    '''
    mode:
        - 0: off
        - 1: closed loop pid
        - 2: zone
        - 3: open loop
        - 4: monitor out
        - 5: warmup supply
    '''

    lsobj, lsobj_passed = get_lsobj(lsobj)
    lsobj.command(f'OUTMODE {output}, {mode}, {INPUTDICT[input]}, {powerup_enable}')
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return None
    else:
        return lsobj

def read_outmode(lsobj=None, output=OUTPUT):

    lsobj, lsobj_passed = get_lsobj(lsobj)
    outmode = lsobj.query(f'OUTMODE? {output}')
    message = [int(i) for i in outmode.split(',')]
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return message, None
    return message, lsobj

def set_mout(mout, lsobj=None, output=OUTPUT, powerup_enable=1):

    lsobj, lsobj_passed = get_lsobj(lsobj)
    lsobj.command(f'MOUT {output}, {mout}')
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return None
    else:
        return lsobj

def read_mout(lsobj=None, output=OUTPUT):

    lsobj, lsobj_passed = get_lsobj(lsobj)
    mout = float(lsobj.query(f'MOUT? {output}'))
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return mout, None
    return mout, lsobj

def set_zone(zone, upper_bound, p, i, d, range, rate, lsobj=None, mout=0, output=OUTPUT, input=INPUT):

    lsobj, lsobj_passed = get_lsobj(lsobj)
    lsobj.command(f'ZONE {output}, {zone}, {upper_bound}, {p}, {i}, {d}, {mout}, {range}, {INPUTDICT[input]}, {rate}')
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return None
    else:
        return lsobj

def read_zone(zone, lsobj=None, output=OUTPUT):

    lsobj, lsobj_passed = get_lsobj(lsobj)
    message = lsobj.query(f'ZONE? {output}, {zone}')
    message = [float(i) for i in message.split(',')]
    if lsobj_passed==False:
        close_lakeshore(lsobj)
        return message, None
    return message, lsobj

def initialize_lakeshore():
    return ls.model_350.Model350()

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
