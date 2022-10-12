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

# confiure system based on config file
with open(os.path.dirname(__file__)+ r'\..\..\Configuration.txt', "r") as f_conf:
    conf_info = f_conf.read()
    conf_info_split = conf_info.split('\n')
    device_id = conf_info_split[0].split('\t')[1]
    port_id = conf_info_split[1].split('\t')[1]

channel_name = ['/%s/demods/0/sample','/%s/demods/1/sample','/%s/demods/2/sample','/%s/demods/3/sample']

def move_x(position, tolerance=1, go_back=10):
    # Attocube initialization
    ax = {'x':0, 'y':1, 'z':2}
    anc = Positioner()

    pos = float(position)
    target = float(pos-go_back)
    tol = float(tolerance)

    # Move to go_back position to prevent hysteresis
    anc.moveAbsolute(ax['x'], int(target*1000))
    time.sleep(0.1)
    error = np.abs(target-anc.getPosition(ax['x'])/1000)
    while (error >= tol):
        anc.moveAbsolute(ax['x'], int(target*1000))
        time.sleep(0.1)
        error = np.abs(target-anc.getPosition(ax['x'])/1000)
    # Move to specified position
    anc.moveAbsolute(ax['x'], int(pos*1000))
    time.sleep(0.1)
    error = np.abs(pos-anc.getPosition(ax['x'])/1000)
    while (error >= tol):
        anc.moveAbsolute(ax['x'], int(pos*1000))
        time.sleep(0.1)
        error = np.abs(pos-anc.getPosition(ax['x'])/1000)

    print(anc.getPosition(ax['x'])/1000)
    anc.close()

def move_y(position, tolerance=1, go_back=10):
    # Attocube initialization
    ax = {'x':0, 'y':1, 'z':2}
    anc = Positioner()

    pos = float(position)
    target = float(pos-go_back)
    tol = float(tolerance)

    # Move to go_back position to prevent hysteresis
    anc.moveAbsolute(ax['y'], int(target*1000))
    time.sleep(0.1)
    error = np.abs(target-anc.getPosition(ax['y'])/1000)
    while (error >= tol):
        anc.moveAbsolute(ax['y'], int(target*1000))
        time.sleep(0.1)
        error = np.abs(target-anc.getPosition(ax['y'])/1000)
    # Move to specified position
    anc.moveAbsolute(ax['y'], int(pos*1000))
    time.sleep(0.1)
    error = np.abs(pos-anc.getPosition(ax['y'])/1000)
    while (error >= tol):
        anc.moveAbsolute(ax['y'], int(pos*1000))
        time.sleep(0.1)
        error = np.abs(pos-anc.getPosition(ax['y'])/1000)

    print(anc.getPosition(ax['y'])/1000)
    anc.close()

def move_z(position, tolerance=1, go_back=10):
    # Attocube initialization
    ax = {'x':0, 'y':1, 'z':2}
    anc = Positioner()

    pos = float(position)
    target = float(pos-go_back)
    tol = float(tolerance)

    # Move to go_back position to prevent hysteresis
    anc.moveAbsolute(ax['z'], int(target*1000))
    time.sleep(0.1)
    error = np.abs(target-anc.getPosition(ax['z'])/1000)
    while (error >= tol):
        anc.moveAbsolute(ax['z'], int(target*1000))
        time.sleep(0.1)
        error = np.abs(target-anc.getPosition(ax['z'])/1000)
    # Move to specified position
    anc.moveAbsolute(ax['z'], int(pos*1000))
    time.sleep(0.1)
    error = np.abs(pos-anc.getPosition(ax['z'])/1000)
    while (error >= tol):
        anc.moveAbsolute(ax['z'], int(pos*1000))
        time.sleep(0.1)
        error = np.abs(pos-anc.getPosition(ax['z'])/1000)

    print(anc.getPosition(ax['z'])/1000)
    anc.close()

def read_position():
    # Attocube initialization
    ax = {'x':0, 'y':1, 'z':2}
    anc = Positioner()

    x=anc.getPosition(ax['x'])/1000
    y=anc.getPosition(ax['y'])/1000
    z=anc.getPosition(ax['z'])/1000

    print('x: '+str(x))
    print('y: '+str(y))
    print('z: '+str(z))

    anc.close()

def rotate_axis(angle, axis_index):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    axis_rot = newport.NewportESP301Axis(controller,axis_index-1)
    axis_rot.enable()

    axis_rot.move(angle,absolute=True)

def read_angle(axis_index):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    axis_rot = newport.NewportESP301Axis(controller,axis_index-1)
    axis_rot.enable()

    print(axis_rot.position)

def set_temperature(temperature):
    '''
    sets lakeshore setpoint

    args:
        - temperature(float)

    returns: None
    '''
    # Lakeshore initialization
    lsobj = ls.initialization_lakeshore335()
    temp=float(temperature)
    ls.set_ramp(lsobj, 1, 0, 0)
    ls.set_setpoint(lsobj, 1, temp)
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

def read_temperature():
    '''
    reads temperature from lakeshore controller

    args: None

    returns:
        - temp(float):  read temperature
    '''
    lsobj = ls.initialization_lakeshore335()
    temp = ls.read_temperature(lsobj)
    ls.close_lakeshore335(lsobj)
    return temp

def z_scan(z_start, z_end, step_size):
    pass

def get_input():
    global button
    while True:
        if button == False:
            break
        else:
            newValue = input('Stop? (Y/N) ')
            time.sleep(1)
            if (newValue == 'Y'):
                button = False
            if (button == False):
                break

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

def save_to_file(fname, data, header, metadata=None):
    '''
    utility function for saving data to a file, with optional metadata

    args:
        - fname(string):           full path to datafile
        - data(array):             (n,m) array containing data
        - header(array):           (m) array of strings labeling each column
        - metadata(str):           in progress

    returns: None
    '''
    if not(len(header) == len(data[0,:])):
        raise ValueError('number of header items does not match number of data columns.')
    with open(fname, 'w') as f:
        for item in header:
            f.write(str(item)+'\t')
        f.write('\n')
        for line in data:
            for item in line:
                f.write(str(item)+'\t')
            f.write('\n')

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
