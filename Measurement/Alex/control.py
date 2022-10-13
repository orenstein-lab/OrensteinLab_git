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
from IPython.display import clear_output
import threading
import ipywidgets as widgets
from pyanc350.v2 import Positioner
import OrensteinLab_git.Instrument.Lakeshore.Lakeshore335 as ls
from strain_control.strain_client import StrainClient
'''
To do:
    - every motor move function takes arguments (setpoint, obj=None, kwargs)
    - for each motor write an initialize function that returns an obj and a close function that closes that motor given an obj.
    - clean up rotate_axis and make rotate_axis_1 and 2 wrapper functions
    - write a set_coil function
    - anything else that seems fitting to go in here!
    - for all of these functions having nice printing feedback would be nice.
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

################
### Medthods ###
################

def move_attocube(axis, position, torlerance=1, go_back=10, anc=None):
    '''
    utility to move attocube
    '''

    # Attocube initialization
    ax = {'x':0, 'y':1, 'z':2}
    anc_passed = True
    if anc==None:
        anc = Positioner()
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
        anc.close()

def read_attocube(anc=None, print=True):
    '''
    Read all attocube positions, printing and returning as desired
    '''
    # Attocube initialization
    ax = {'x':0, 'y':1, 'z':2}
    anc_passed = True
    if anc==None:
        anc = Positioner()
        anc_passed = False

    time.sleep(0.1)
    x=anc.getPosition(ax['x'])/1000
    y=anc.getPosition(ax['y'])/1000
    z=anc.getPosition(ax['z'])/1000

    if print==True:
        print('x: '+str(x))
        print('y: '+str(y))
        print('z: '+str(z))

    # print and close only if another process hasn't passed anc object
    if anc_passed == False:
        print(anc.getPosition(ax[axis])/1000)
        anc.close()

    return x, y, z

def move_x(position, tolerance=1, go_back=10, anc=None):
    '''
    wrapper to move attocube x positioner.
    '''
    move_attocube('x', position, tolerance, go_back, anc)

def move_y(position, tolerance=1, go_back=10, anc==None):
    '''
    wrapper to move attocube y positioner.
    '''
    move_attocube('y', position, tolerance, go_back, anc)

def move_z(position, tolerance=1, go_back=10):
    '''
    wrapper to move attocube z positioner.
    '''
    move_attocube('z', position, tolerance, go_back, anc)

def rotate_axis(axis_index, angle, controller=None, axis_rot=None):
    '''
    rotate waveplate. I can now add any checks that have to happen typically.
    '''
    #ESP301 initialization
    if controller==None:
        controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # initialize axis
    if axis_rot==None:
        axis_rot = newport.NewportESP301Axis(controller,axis_index-1)
        axis_rot.enable()

    axis_rot.move(angle,absolute=True)
    check_axis_stability(axis_rot)

def check_axis_stability(axis_rot):
    '''
    helper function for rotate_axis.
    '''
    while True:
        time.sleep(0.03)
        try:
            if axis_rot.is_motion_done==True:
                break
        except ValueError:
            pass

def read_angle(axis_index, controller=None, axis_rot=None, print=True):
    '''
    read angle on an axis
    '''
    # ESP301 initialization
    if controller==None:
        controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # initialize axis
    if axis_rot==None:
        axis_rot = newport.NewportESP301Axis(controller,axis_index-1)
        axis_rot.enable()

    pos = axis_rot.position
    if print==True:
        print(pos)
    return pos

def set_temperature(temperature, tolerance=0.01, wait_time=0, max_check=0, lsobj=None):
    '''
    sets lakeshore setpoint, waits until temperature is within tolerance of setpoint, and waits for soak time before returning.

    args:
        - temperature(float)

    returns: None
    '''
    # Lakeshore initialization
    lsobj_passed=True
    if lsobj==None:
        lsobj = ls.initialization_lakeshore335()
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
        ls.close_lakeshore335(lsobj)

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

def get_temp_values(starttemp, endtemp, stepsize=0.):
    '''Returns an array of temperature values for a given starting temperature (starttemp), ending temperature (endtemp),
    and step size between temperature points (stepsize).
    starttemp: value to start list of temperatures at (in units of K)
    endtemp: maximum value for temperature to ramp up to (in units of K)
    stepsize: amount to increment temp by during ramp (in units of K)'''

    temprange = endtemp-starttemp
    if stepsize == 0:
        return np.asarray(starttemp)
    else:
        roundsteps = np.rint([temprange/stepsize])[0]
    if temprange/stepsize != roundsteps:
        newendtemp = starttemp + roundsteps*stepsize
        print('Warning: Specified end temperature value cannot be reached given the specified start temperature and step size. '
              +'End temperature will be rounded to the nearest step. New end temperature: '+str(newendtemp)+'K.')
        return np.linspace(starttemp, newendtemp, int(roundsteps)+1)
    else:
        return np.linspace(starttemp, endtemp, int(roundsteps)+1)

def get_angle_range(startangle, endangle, stepsize=0.):
    '''Returns an array of angle values for a given starting angle (startangle), ending angle (endangle),
    and step size between angle points (stepsize).
    startangle: value to start list of angles at (in units of deg)
    endangle: maximum value for angle to ramp up to (in units of deg)
    stepsize: amount to increment angle by during ramp (in units of deg)'''

    anglerange = endangle-startangle
    if stepsize == 0:
        return np.asarray(startangle)
    else:
        roundsteps = np.rint([anglerange/stepsize])[0]
    if anglerange/stepsize != roundsteps:
        newendangle = startangle + roundsteps*stepsize
        print('Warning: Specified end angle value cannot be reached given the specified start angle and step size. '
              +'End angle will be rounded to the nearest step. New end angle: '+str(newendangle)+' deg.')
        return np.linspace(startangle, newendangle, int(roundsteps)+1)
    else:
        return np.linspace(startangle, endangle, int(roundsteps)+1)

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
