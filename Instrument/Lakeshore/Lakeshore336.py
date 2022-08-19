#!/usr/bin/env python
# coding: utf-8
# @author: Yue Sun, UCB


import lakeshore

def initialization_lakeshore336():
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

def close_lakeshore336(inst):
    inst.disconnect_usb()

