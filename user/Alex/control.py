'''
Main file for controlling lab equipment and orchestrating measurements
'''

get_ipython().run_line_magic('matplotlib', 'notebook')
import zhinst.utils as ziutils
import zhinst.core as zicore
import instruments.newport as newport
from pymeasure.instruments.newport import ESP300
import OrensteinLab_git.Instrument.OptiCool.OptiCool_Control as optc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os.path
import time
import datetime
import math
import pickle
from IPython.display import clear_output
import threading
import ipywidgets as widgets
from pyanc350.v2 import Positioner
import OrensteinLab_git.Instrument.Lakeshore.Lakeshore as ls
from strain_control.strain_client import StrainClient





def set_opticool_temperature(temperature, rate=10, mode=0, obj=None):

    obj_passed=True
    if obj==None:
        obj = initialize_opticool()
        obj_passed==False

    optc.set_temperature(obj, temperature, rate, mode)

    if obj_passed==False:
        close_opticool(obj)

def read_opticool_temperature(obj=None):

    obj_passed=True
    if obj==None:
        obj = initialize_opticool()
        obj_passed==False

    temp = optc.read_temperature(obj)[0]

    if obj_passed==False:
        close_opticool(obj)

    return temp

def set_opticool_field(field, rate=110, mode=0, obj=None):

    obj_passed=True
    if obj==None:
        obj = initialize_opticool()
        obj_passed==False

    optc.set_field(obj, field, rate, mode)

    if obj_passed==False:
        close_opticool(obj)

def read_opticool_field(obj=None):

    obj_passed=True
    if obj==None:
        obj = initialize_opticool()
        obj_passed==False

    field = optc.read_field(obj)[0]

    if obj_passed==False:
        close_opticool(obj)

    return field

#######################################
### Wrapper Move and Read Functions ###
#######################################

def initialize_opticool():
    return optc.connect_opticool()

#####################
### Close Methods ###
#####################


def close_opticool(obj):
    optc.disconnect_opticool(obj)
