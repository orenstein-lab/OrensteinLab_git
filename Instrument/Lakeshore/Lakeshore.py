#!/usr/bin/env python
# coding: utf-8
# @author: Yue Sun, UCB


import lakeshore

def initialize_lakeshore335():
    return lakeshore.model_335.Model335(57600)

def initialize_lakeshore336():
    return lakeshore.model_336.Model336()

def read_temperature(inst):
    return float(inst.query('KRDG?'))

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

def close_lakeshore(inst):
    inst.disconnect_usb()
