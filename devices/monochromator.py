#!/usr/bin/env python
# coding: utf-8
# @author: Yue Sun, UCB

# from OrensteinLab_git.configuration import config_dict
import time
import numpy as np
import sys

try:
  import pyvisa
except ImportError:
  # Change this to match your pyvisa installation path.
  sys.path.append( r"C:\Users\BigLab2020\anaconda3\lib\site-packages\pyvisa" )
  import pyvisa

def read_wavelength(obj=None):

    obj, obj_passed = get_obj(obj)
    wavelength_ls = obj.query('WAVE?').split()
    wavelength = float(wavelength_ls[0])
    if obj_passed==False:
        print(wavelength)
        close_monochromator(obj)
        return wavelength, None
    else:
        return wavelength, obj

def read_filter(obj=None):

    obj, obj_passed = get_obj(obj)
    filter_ls = obj.query('FILTER?').split()
    filter = int(filter_ls[0])
    if obj_passed==False:
        print(filter)
        close_monochromator(obj)
        return filter, None
    else:
        return filter, obj

def read_grating(obj=None):

    obj, obj_passed = get_obj(obj)
    grating_ls = obj.query('GRATing?')
    grating = int(grating_ls[0][0])
    if obj_passed==False:
        print(grating)
        close_monochromator(obj)
        return grating, None
    else:
        return grating, obj

def set_filter(filter_ind, obj=None):

    obj, obj_passed = get_obj(obj)
    command = 'FILTER '+str(filter_ind)
    obj.write(command)
    if obj_passed==False:
        close_monochromator(obj)
        return None
    else:
        return obj

def set_grating(grating_ind, obj=None):

    obj, obj_passed = get_obj(obj)
    command = 'GRATing '+str(grating_ind)
    obj.write(command)
    if obj_passed==False:
        close_monochromator(obj)
        return None
    else:
        return obj

def set_wavelength(wavelength, obj=None):

    obj, obj_passed = get_obj(obj)
    filter_ind = 0

    if wavelength < 800:
        filter_ind = 1 # empty slot: 450 ~ 800 nm
    elif wavelength < 1200:
        filter_ind = 2 # 10CGA-780 (LP filter), 800 ~ 1200 nm
    elif wavelength < 1900:
        filter_ind = 3 # 10CGA-1000 (LP filter), 1200 ~ 1900 nm
#     elif wavelength < 2200: # 90107911, 1900 ~ 2200 nm
#         filter_ind = 4
    else:
        ValueError('Wavelength out of range!')

    if filter_ind > 0:
        filter_ind_old = read_filter(obj)
        if filter_ind != filter_ind_old :
            set_filter(filter_ind,obj)
            time.sleep(1)
            filter_ind_new = read_filter(obj)
            if filter_ind_new != filter_ind:
                ValueError('Fail to change filter')

    # the following grating selection is for CS130B-3-MC
#     if wavelength < 650:
#         set_grating(1,obj)
#     else:
#         set_grating(2,obj)
    command = 'GOWAVE '+str(round(wavelength,3))
    obj.write(command)
    if obj_passed==False:
        close_monochromator(obj)
        return None
    else:
        return obj

def initialize_monochromator():
    found_unit = False
    rm = pyvisa.ResourceManager()
    instr_list = rm.list_resources()
    for instr in instr_list :
        fields = instr.split('::')
        if( len(fields) == 5 ):         # 5 fields in a usbtmc device entry of the instrument list.
            if( fields[0] == "USB0" ):      # ... should probably not care about the '0'...
                if( fields[1] == "0x1FDE" ):    # Newport
                    if(( fields[2] == "0x0014" ) or( fields[2] == "0x0006")):    # CS260B or CS130B
                        addr_str = instr
                        found_unit = True
                        break   # Found one! Stop examining the instrument list

    if found_unit:
    #    print(addr_str)
       return rm.open_resource(addr_str)
    else:
       raise ValueError('Cannot Find Monochromator')

def close_monochromator(obj):
    obj.close()

def get_obj(obj):
    '''
    helper function for getting object both within another script or on the command line.
    '''
    obj_passed = True
    if obj==None:
        obj = initialize_monochromator()
        obj_passed=False
    return obj, obj_passed
