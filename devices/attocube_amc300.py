import numpy as np
from pyamc300 import AMC
from OrensteinLab_git.configuration import config_dict
import pickle
import dill
import time
import os

axis_dict = {'x':0, 'y':1, 'z':2}
ATTOCUBE_HANDLE_FNAME = config_dict['Attocube Handle']
ATTOCUBE_IP = config_dict['AMC300 IP']

######################
### Core Functions ###
######################

def move_attocube(axis, position, amc=None, tolerance=1, go_back=0, ground=False):
    '''
    utility to move attocube
    '''

    # Attocube initialization
    #amc_passed = True
    # if amc==None:
    amc = initialize_attocube()
    amc_passed = False
    go_back=10


    pos = float(position)
    target = float(pos-go_back)
    tol = float(tolerance)

    # Move to go_back position to prevent hysteresis
    # set_attocube_output(axis,True)
    amc.control.setControlOutput(axis_dict[axis], True)
    amc.move.setControlTargetPosition(axis_dict[axis], target * 1e3)
    amc.control.setControlMove(axis_dict[axis], True)

    while not amc.status.getStatusTargetRange(axis):
        # # Read out position in nm
        if not amc.status.getStatusTargetRange(axis):
            position_real = amc.move.getPosition(axis)
            amc.control.setControlMove(axis, True)
            time.sleep(0.1)
        if repeat_time > 10:
            amc.move.setControlContinuousBkwd(axis, True)
            time.sleep(0.2)
            #    # Stop
            amc.move.setControlContinuousBkwd(axis, False)
            repeat_time = 0

    # Move to specified position
    amc.move.setControlTargetPosition(axis, pos * 1e3)
    amc.control.setControlMove(axis, True)
    while not amc.status.getStatusTargetRange(axis):
        # # Read out position in nm
        if not amc.status.getStatusTargetRange(axis):
            position_real = amc.move.getPosition(axis)
            amc.control.setControlMove(axis, True)
            time.sleep(0.1)
        if repeat_time > 10:
            amc.move.setControlContinuousBkwd(axis, True)
            time.sleep(0.2)
            #    # Stop
            amc.move.setControlContinuousBkwd(axis, False)
            repeat_time = 0        
    '''
    amc.moveAbsolute(axis_dict[axis], int(target/1e6))
    time.sleep(0.1)
    error = np.abs(target-anc.getPosition(axis_dict[axis])*1e6)
    while (error >= tol):
        anc.moveAbsolute(axis_dict[axis], int(target/1e6))
        time.sleep(0.1)
        error = np.abs(target-anc.getPosition(axis_dict[axis])*1e6)

    # Move to specified position
    anc.moveAbsolute(axis_dict[axis], int(pos/1e6))
    time.sleep(0.1)
    error = np.abs(pos-anc.getPosition(axis_dict[axis])*1e6)
    while (error >= tol):
        anc.moveAbsolute(axis_dict[axis], int(pos/1e6))
        time.sleep(0.1)
        error = np.abs(pos-anc.getPosition(axis_dict[axis])*1e6)
    '''
    # print and close only if another process hasn't passed anc object
    if amc_passed == False:
        print(amc.move.getPosition(axis_dict[axis])*1e-3)
        close_attocube(amc)

def read_attocube(axis, amc=None, print_flag=True):
    '''
    Read all attocube positions, printing and returning as desired
    '''
    # Attocube initialization
    amc_passed = True
    if amc==None:
        amc = initialize_attocube()
        amc_passed = False

    if axis not in list(axis_dict.keys()):
        raise ValueError('Invalid axis, please choose from ["x", "y", "z"].')

    time.sleep(0.1)
    pos=amc.move.getPosition(axis_dict[axis])*1e-3

    if print_flag==True and amc_passed==False:
        print(f'{axis}: {pos}')

    # print and close only if another process hasn't passed anc object
    if amc_passed == False:
        print(pos)
        close_attocube(amc)

    return pos

def set_attocube_output(axis, state, amc=None):
    '''
    ground or unground attocube
    '''

    amc_passed = True
    if amc==None:
        amc = initialize_attocube()
        amc_passed = False
    
    amc.control.setControlOutput(axis_dict[axis], state)

    if amc_passed == False:
        print(amc.move.getPosition(axis_dict[axis])*1e-3)
        close_attocube(amc)

def initialize_attocube():
    '''
    attocube initialization function, which checks if an initialization already exists at attocube_handle before initalizizing another instance.
    '''

    amc = AMC.Device(ATTOCUBE_IP)
    amc.connect()
    '''
    try:
        if os.path.exists(ATTOCUBE_HANDLE_FNAME):
            amc = AMC.Device(ATTOCUBE_IP)
            amc.connect()
            print('death')
           # with open(ATTOCUBE_HANDLE_FNAME, 'rb') as f:
                # amc = dill.load(f)
                #read_attocube('x',amc) # test the connection
        else:
            amc = AMC.Device(ATTOCUBE_IP)
            amc.connect()
            #with open(ATTOCUBE_HANDLE_FNAME, 'wb') as f:
                #dill.dump(amc, f)
    except Exception as e:
        print(f'initialization error: {e}')
        os.remove(ATTOCUBE_HANDLE_FNAME)
        amc = AMC.Device(ATTOCUBE_IP)
        amc.connect()
        #with open(ATTOCUBE_HANDLE_FNAME, 'wb') as f:
            #dill.dump(amc, f)
    '''
    return amc

def close_attocube(amc):
        if os.path.exists(ATTOCUBE_HANDLE_FNAME):
            os.remove(ATTOCUBE_HANDLE_FNAME)
            amc.close()

#########################
### Wrapper Functions ###
#########################

def move_pos(x, y, z=0):
    move_x(x)
    move_y(y)
    if z==0:
        pass
    else:
        move_z(z)

def move_x(position,  amc=None, tolerance=1, go_back=0):
    '''
    wrapper to move attocube x positioner.
    '''
    move_attocube('x', position, amc, tolerance, go_back)

def move_y(position,  amc=None, tolerance=1, go_back=0):
    '''
    wrapper to move attocube y positioner.
    '''
    move_attocube('y', position, amc, tolerance, go_back)

def move_z(position,  amc=None, tolerance=1, go_back=0):
    '''
    wrapper to move attocube z positioner.
    '''
    move_attocube('z', position, amc, tolerance, go_back)

def read_x(amc=None, print_flag=True):
    return read_attocube('x', amc, print_flag)

def read_y(amc=None, print_flag=True):
    return read_attocube('y', amc, print_flag)

def read_z(amc=None, print_flag=True):
    return read_attocube('z', amc, print_flag)
