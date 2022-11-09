'''
Main file for controlling lab equipment and orchestrating measurements
'''

get_ipython().run_line_magic('matplotlib', 'notebook')
import zhinst.utils as ziutils
import instruments.newport as newport
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
import OrensteinLab_git.Instrument.Lakeshore.Lakeshore335 as ls
from strain_control.strain_client import StrainClient
'''
To do:
    - clean up printing to eliminate erroneous prints/only print when called from the command line or when we intuitively want it to print
    - every motor move function takes arguments (setpoint, obj=None, kwargs)
    - write a set_coil function
    - anything else that seems fitting to go in here!
'''

#####################
### Configuration ###
#####################
with open(os.path.dirname(__file__)+ r'\..\..\Configuration.txt', "r") as f_conf:
    conf_info = f_conf.read()
    conf_info_split = conf_info.split('\n')
    device_id = conf_info_split[0].split('\t')[1]
    port_id = conf_info_split[1].split('\t')[1]
channel_name = ['/%s/demods/0/sample','/%s/demods/1/sample','/%s/demods/2/sample','/%s/demods/3/sample']
ATTOCUBE_HANDLE_FNAME = f'{os.path.dirname(__file__)}\\..\\..\\attocube_handle'

###############################
### Instrument Read Methods ###
###############################

def read_zurich_lockin(daq_objs=None, time_constant=0.3, channel_index=1, R_channel_index=1):

    # initialize
    if daq_objs==None:
        daq, device, props = ctrl.initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    time.sleep(time_constant*4)
    sample = daq.getSample(channel_name[channel_index-1] % device)
    sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
    x = sample["x"][0]
    y = sample["y"][0]
    r = sample["R"][0]
    sample_R = daq.getSample(channel_name[R_channel_index - 1] % device)
    sample_R["R"] = np.abs(sample_R["x"] + 1j * sample_R["y"])
    x_R = sample_R["x"][0]
    y_R = sample_R["y"][0]
    r_R = sample_R["R"][0]

    return x, y, r, x_R, y_R, r_R

###################################
### Motor Move and Read Methods ###
###################################

def move_attocube(axis, position, anc=None, tolerance=1, go_back=10):
    '''
    utility to move attocube
    '''

    # Attocube initialization
    ax = {'x':0, 'y':1, 'z':2}
    anc_passed = True
    if anc==None:
        anc = initialize_attocube()
        anc_passed = False


    pos = float(position)
    target = float(pos-go_back)
    tol = float(tolerance)

    # Move to go_back position to prevent hysteresis
    anc.moveAbsolute(ax[axis], int(target*1000))
    time.sleep(0.1)
    error = np.abs(target-anc.getPosition(ax[axis])/1000)
    while (error >= tol):
        anc.moveAbsolute(ax[axis], int(target*1000))
        time.sleep(0.1)
        error = np.abs(target-anc.getPosition(ax[axis])/1000)

    # Move to specified position
    anc.moveAbsolute(ax[axis], int(pos*1000))
    time.sleep(0.1)
    error = np.abs(pos-anc.getPosition(ax[axis])/1000)
    while (error >= tol):
        anc.moveAbsolute(ax[axis], int(pos*1000))
        time.sleep(0.1)
        error = np.abs(pos-anc.getPosition(ax[axis])/1000)

    # print and close only if another process hasn't passed anc object
    if anc_passed == False:
        print(anc.getPosition(ax[axis])/1000)
        close_attocube(anc)

def read_attocube(axis, anc=None, print_flag=True):
    '''
    Read all attocube positions, printing and returning as desired
    '''
    # Attocube initialization
    ax_dict = {'x':0, 'y':1, 'z':2}
    anc_passed = True
    if anc==None:
        anc = initialize_attocube()
        anc_passed = False

    if axis not in list(ax_dict.keys()):
        raise ValueError('Invalid axis, please choose from ["x", "y", "z"].')

    time.sleep(0.1)
    pos=anc.getPosition(ax_dict[axis])/1000

    if print_flag==True and anc_passed==False:
        print(f'{axis}: {pos}')

    # print and close only if another process hasn't passed anc object
    if anc_passed == False:
        print(pos)
        close_attocube(anc)

    return pos

def move_x(position,  anc=None, tolerance=1, go_back=10):
    '''
    wrapper to move attocube x positioner.
    '''
    move_attocube('x', position, anc, tolerance, go_back)

def move_y(position,  anc=None, tolerance=1, go_back=10):
    '''
    wrapper to move attocube y positioner.
    '''
    move_attocube('y', position, anc, tolerance, go_back)

def move_z(position,  anc=None, tolerance=1, go_back=10):
    '''
    wrapper to move attocube z positioner.
    '''
    move_attocube('z', position, anc, tolerance, go_back)

def read_x(anc=None, print_flag=True):
    return read_attocube('x', anc, print_flag)

def read_y(anc=None, print_flag=True):
    return read_attocube('y', anc, print_flag)

def read_z(anc=None, print_flag=True):
    return read_attocube('z', anc, print_flag)

def rotate_axis(axis_index, angle, axis):
    '''
    rotate waveplate. I can now add any checks that have to happen typically.
    '''
    # initialize axis
    if axis==None:
        axis = initialize_rotation_axis(axis_index)

    axis.move(angle,absolute=True)
    check_axis_stability(axis)

def read_axis(axis_index, axis=None, print_flag=True):
    '''
    read angle on an axis
    '''
    # initialize axis
    obj_passed = True
    if axis==None:
        axis = initialize_rotation_axis(axis_index)
        obj_passed = False

    while True:
        time.sleep(0.03)
        try:
            pos = float(axis.position)
            break
        except:
            pass

    if print_flag==True and obj_passed==False:
        print(pos)
    return pos

def check_axis_stability(axis):
    '''
    helper function for rotate_axis.
    '''
    while True:
        time.sleep(0.03)
        try:
            if axis.is_motion_done==True:
                break
        except:
            pass

def rotate_axis_1(angle, axis=None):
    rotate_axis(1, angle, axis)

def rotate_axis_2(angle, axis=None):
    rotate_axis(2, angle, axis=None)

def read_axis_1(axis=None, print_flag=True):
    return read_axis(1, axis, print_flag)

def read_axis_2(axis=None, print_flag=True):
    return read_axis(2, axis, print_flag)

def set_temperature(temperature, lsobj=None, tolerance=0.01, wait_time=0, max_check=300):
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
    ls.set_ramp(lsobj, 1, 0, 0)
    ls.set_setpoint(lsobj, 1, temp)
    time.sleep(0.1)

    # check for stability
    current_temp = []
    for m in range(max_check):
        current_temp.append(ls.read_temperature(lsobj))
        if m >= 2 and abs(np.mean(current_temp[-3:]) - temp) < 0.05:
            time.sleep(wait_time)
            break
        else:
            time.sleep(1)

    if lsobj_passed==False:
        close_lakeshore(lsobj)

def set_lakeshore_range(range):
    '''
    sets lakeshore range (0=off, 1=Low, 2=Med, 3=High)

    args:
        - range(int)

    returns: None
    '''
    range=int(range)
    if range not in [0,1,2,3]:
        raise ValueError(f'{range} is not a valid range. Please choose from (0=off, 1=Low, 2=Med, 3=High)')
    lsobj = ls.initialization_lakeshore335()
    ls.set_range(lsobj, 1, range)
    ls.close_lakeshore335(lsobj)

def read_temperature(lsobj=None):
    '''
    reads temperature from lakeshore controller

    args: None

    returns:
        - temp(float):  read temperature
    '''
    lsobj_passed = True
    if lsobj==None:
        lsobj = ls.initialization_lakeshore335()
        lsobj_passed=False
    temp = ls.read_temperature(lsobj)
    if lsobj_passed==False:
        ls.close_lakeshore335(lsobj)
    return temp

def set_coil():
    return 1

def read_coil():
    return 1

def set_strain_voltage(channel, voltage, sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    sc.set_voltage(channel, voltage)

def read_strain_capacitance(sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    cap = sc.get_cap()

    if sc_passed == False:
        print(pos)

    return cap

def set_strain_voltage_compress(voltage, sc=None):
    set_strain_voltage(2,voltage, sc=sc)

def set_strain_voltage_tension(voltage, sc=None):
    set_strain_voltage(2,voltage, sc=sc)

###############################
### Initialization Medthods ###
###############################

def initialize_zurich_lockin():
    apilevel=6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)
    return daq, device, props

def initialize_esp():
    return newport.NewportESP301.open_serial(port=port_id, baud=921600)

def initialize_attocube():
    # can this function check if the attocube has already been initialized?
    try: # if there is already a connection, this fails and we search for the handle in file.
        anc = Positioner()
        with open(ATTOCUBE_HANDLE_FNAME, 'wb') as f:
            pickle.dump(anc, f)
    except:
        with open(ATTOCUBE_HANDLE_FNAME, 'rb') as f:
            anc = pickle.load(f)
    return anc

def initialize_rotation_axis(index):
    controller = initialize_esp()
    while True:
        try:
            axis_rot = newport.NewportESP301Axis(controller, index-1)
            axis_rot.enable()
            break
        except:
            time.sleep(0.1)
    return axis_rot

def initialize_rot_axis_1():
    return initialize_rotation_axis(1)

def initialize_rot_axis_2():
    return initialize_rotation_axis(2)

def initialize_lakeshore():
    return ls.initialization_lakeshore335()

def initialize_strain_cell_client():
    return StrainClient()

#####################
### Close Methods ###
#####################

def close_attocube(anc):
    try:
        anc.close()
    except:
        return 0

def close_lakeshore(obj):
    ls.close_lakeshore335(obj)

def close_rot_axis_1(obj):
    return 0

def close_rot_axis_2(obj):
    return 0

def close_lockin(obj):
    return 0

def close_strain_cell_client(obj):
    return 0

''' a script I wrote for collecting data on capacitive sensor
filename = 'G:/Shared drives/Orenstein Lab/Data/alex/20220818_strain_cell_capacitor_noise_25ft_coax_300kHz.txt'
num_points = 1000
t = np.zeros(num_points)
cap = np.zeros(num_points)
t0 = time.time()
i = 0
t_old = t0
while i < num_points:
    t_new = time.time()
    if t_new - t_old > 0.1:
        c = sc.get_cap()
        cap[i] = c
        t[i] = t_new - t0
        i = i + 1
        t_old = t_new
save_to_file(filename, np.transpose([t, cap]), ['Time', 'Cap'])
'''
