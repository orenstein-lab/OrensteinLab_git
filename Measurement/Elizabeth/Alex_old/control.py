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
import OrensteinLab_git.Instrument.Lakeshore.Lakeshore335 as ls
from strain_control.strain_client import StrainClient
'''
To do:
    - clean up printing to eliminate erroneous prints/only print when called from the command line or when we intuitively want it to print
    - every motor move function takes arguments (setpoint, obj=None, kwargs)
    - write a set_coil function
    - anything else that seems fitting to go in here!

to find pymeasure instruments:

from pymeasure.instruments import list_resources
list_resources()
'''

#####################
### Configuration ###
#####################
with open(os.path.dirname(__file__)+ r'\..\..\Configuration.txt', "r") as f_conf:
    conf_info = f_conf.read()
    conf_info_split = conf_info.split('\n')
    device_id = conf_info_split[0].split('\t')[1]
    port_id = conf_info_split[1].split('\t')[1]
    esp300_id = conf_info_split[2].split('\t')[1]
ATTOCUBE_HANDLE_FNAME = f'{os.path.dirname(__file__)}\\..\\..\\attocube_handle'

###############################
### Instrument Read Methods ###
###############################

def read_zurich_lockin(daq_objs=None, time_constant=0.3, poll_timeout=500, channel_index=1, R_channel_index=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/demods/{channel_index-1}/timeconstant' % device, time_constant)
    daq.setDouble(f'/%s/demods/{R_channel_index-1}/timeconstant' % device, time_constant)
    poll_length = time_constant
    poll_timeout = poll_timeout # ms
    poll_flags = 0
    poll_return_flat_dict = True

    # subscribe to channels and read mfli
    time.sleep(time_constant*4)
    daq.subscribe(f'/%s/demods/{channel_index-1}/sample' % device)
    daq.subscribe(f'/%s/demods/{R_channel_index-1}/sample' % device)
    mfli_dict = daq.poll(poll_length, poll_timeout, poll_flags, poll_return_flat_dict)
    data_dict = mfli_dict[f'/%s/demods/{channel_index-1}/sample' % device]
    R_dict = mfli_dict[f'/%s/demods/{R_channel_index-1}/sample' % device]
    daq.unsubscribe('*')

    # extract data from mfli_dict
    x = data_dict['x']
    y = data_dict['y']
    r = np.abs(x + 1j*y)
    x_avg = np.mean(x)
    y_avg = np.mean(y)
    r_avg = np.mean(r)
    x_R = R_dict['x']
    y_R = R_dict['y']
    r_R = np.abs(x_R + 1j*y_R)
    x_R_avg = np.mean(x_R)
    y_R_avg = np.mean(y_R)
    r_R_avg = np.mean(r_R)

    return x_avg, y_avg, r_avg, x_R_avg, y_R_avg, r_R_avg

########################################
### Core Motor Move and Read Methods ###
########################################

def move_attocube(axis, position, anc=None, tolerance=1, go_back=10):
    '''
    utility to move attocube
    '''

    # Attocube initialization
    ax = {'x':0, 'y':1, 'z':2}
    anc_passed = True
    if anc!=None:
        go_back=0
    if anc==None:
        anc = initialize_attocube()
        anc_passed = False

    pos = float(position)
    target = float(pos-go_back)
    tol = float(tolerance)

    # Move to go_back position to prevent hysteresis if go_back!=0
    if go_back!=0:
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

def move_pos(x, y, z=0):
    move_x(x)
    move_y(y)
    if z==0:
        pass
    else:
        move_z(z)

def rotate_axis(axis_index, angle, axis=None):
    '''
    rotate waveplate. I can now add any checks that have to happen typically.
    '''
    # initialize axis
    if axis==None:
        axis = initialize_rotation_axis(axis_index)

    axis.move(angle,absolute=True)
    check_axis_stability(axis)

def corotate_axes(axis_1_index, axis_2_index, angle_1, angle_2, axis_1=None, axis_2=None):

    if axis_1==None:
        axis_1 = initialize_rotation_axis(axis_1_index)
    if axis_2==None:
        axis_2 = initialize_rotation_axis(axis_2_index)

    axis_1.move(angle_1, absolute=True)
    axis_2.move(angle_2, absolute=True)
    while True:
        time.sleep(0.5)
        try:
            if axis_1.is_motion_done==True and axis_2.is_motion_done==True:
                break
        except:
            print('failed to check axis stability, trying agian.')
            #pass

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
        time.sleep(0.1)
        try:
            pos = float(axis.position)
            break
        except:
            print('failed to read axis, trying again.')
            #pass

    if print_flag==True and obj_passed==False:
        print(pos)
    return pos

def check_axis_stability(axis):
    '''
    helper function for rotate_axis.
    '''
    while True:
        time.sleep(0.1)
        try:
            if axis.is_motion_done==True:
                break
        except:
            print('failed to check axis stability, trying agian.')
            #pass

def move_esp300(axis_index, position, axis=None):

    if axis==None:
        axis = initialize_esp300(axis_index)

    axis.position = position
    while True:
        time.sleep(0.1)
        try:
            if axis.motion_done==True:
                break
        except:
            print('failed to check axis stability, trying agian.')
            #pass

def read_esp300(axis_index, axis=None, print_flag=True):
    obj_passed = True
    if axis==None:
        axis = initialize_esp300(axis_index)
        obj_passed = False

    while True:
        time.sleep(0.1)
        try:
            pos = float(axis.position)
            break
        except:
            print('failed to read axis, trying again.')
            #pass

    if print_flag==True and obj_passed==False:
        print(pos)
    return pos

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
    ls.set_ramp(lsobj, 1, 0, 0)
    ls.set_setpoint(lsobj, 1, temp)
    time.sleep(0.1)

    # check for stability
    current_temp = []
    for m in range(max_check):
        current_temp.append(ls.read_temperature(lsobj))
        if m >= 10 and abs(np.mean(current_temp[-3:]) - temp) < tolerance:
            time.sleep(wait_time)
            break
        else:
            time.sleep(1)
    if m==max_check-1:
        print(f'Maximum time exceeded. Temperature: {ls.read_temperature(lsobj)}')

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

def set_strain_ps(voltage, sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    sc.set_ps(voltage)

def set_strain_capacitance(cap, sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    sc.set_cap(cap)

def read_strain_ps(sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    voltage = sc.get_ps()

    if sc_passed == False:
        print(pos)

    return voltage

def read_strain_capacitance(sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    cap = sc.get_cap()

    if sc_passed == False:
        print(pos)

    return cap

def set_zurich_output_amplitude(amplitude, daq_objs=None, wait_time=0, channel=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/sigouts/0/amplitudes/{channel-1}'%device, amplitude)
    time.sleep(wait_time)

def read_zurich_output_amplitude(daq_objs=None, channel=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    amplitude = daq.getDouble(f'/%s/sigouts/0/amplitudes/{channel-1}'%device)
    return amplitude

def set_zurich_frequency(freq, daq_objs=None, wait_time=0, osc=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/oscs/{osc-1}/freq'%device, freq)
    time.sleep(wait_time)

def read_zurich_frequency(daq_objs=None, osc=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    freq = daq.getDouble(f'/%s/oscs/{osc-1}/freq'%device)
    return freq

#######################################
### Wrapper Move and Read Functions ###
#######################################

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

def rotate_axis_1(angle, axis=None):
    rotate_axis(1, angle, axis)

def rotate_axis_2(angle, axis=None):
    rotate_axis(2, angle, axis)

def read_axis_1(axis=None, print_flag=True):
    return read_axis(1, axis, print_flag)

def read_axis_2(axis=None, print_flag=True):
    return read_axis(2, axis, print_flag)

def move_delay_stage(position, axis=None):
    move_esp300(1, position, axis)

def read_delay_stage(axis=None, print_flag=True):
    return read_esp300(1, axis, print_flag)

###############################
### Initialization Medthods ###
###############################

def initialize_zurich_lockin():
    apilevel=6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)
    return daq, device, props

def initialize_esp():
    return newport.NewportESP301.open_serial(port=port_id, baud=921600)

def initialize_esp300(axis_index):
    obj = ESP300(esp300_id)
    if axis_index==1:
        axis = obj.x
    elif axis_index==2:
        axis = obj.y
    elif axis_index==3:
        axis = obj.phi
    else:
        ValueError('could not initialize ESP300. Please choose axis 1, 2, or 3')
    axis.enable()
    return axis

def initialize_delay_stage():
    return initialize_esp300(1)

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

def close_delay_stage(obj):
    return 0

def close_lockin(obj):
    return 0

def close_strain_cell_client(obj):
    return 0

def close_zurich_lockin(obj):
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
