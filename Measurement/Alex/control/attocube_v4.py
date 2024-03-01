import numpy as np
from pyanc350.v4 import Positioner
from OrensteinLab_git.configuration import config_dict
import pickle
import time
import os

axis_dict = {'x':0, 'y':1, 'z':2}
ATTOCUBE_HANDLE_FNAME = config_dict['Attocube Handle']

######################
### Core Functions ###
######################

def move_attocube(axis, position, anc=None, tolerance=1, go_back=0, ground=False):
    '''
    utility to move attocube
    '''

    # Attocube initialization
    anc_passed = True
    if anc==None:
        anc = initialize_attocube()
        anc_passed = False
        go_back=10


    pos = float(position)
    target = float(pos-go_back)
    tol = float(tolerance)

    # Move to go_back position to prevent hysteresis
    anc.moveAbsolute(axis_dict[axis], int(target/1e6))
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

    # print and close only if another process hasn't passed anc object
    if anc_passed == False:
        print(anc.getPosition(axis_dict[axis])*1e6)
        close_attocube(anc)

def read_attocube(axis, anc=None, print_flag=True):
    '''
    Read all attocube positions, printing and returning as desired
    '''
    # Attocube initialization
    anc_passed = True
    if anc==None:
        anc = initialize_attocube()
        anc_passed = False

    if axis not in list(axis_dict.keys()):
        raise ValueError('Invalid axis, please choose from ["x", "y", "z"].')

    time.sleep(0.1)
    pos=anc.getPosition(axis_dict[axis])*1e6

    if print_flag==True and anc_passed==False:
        print(f'{axis}: {pos}')

    # print and close only if another process hasn't passed anc object
    if anc_passed == False:
        print(pos)
        close_attocube(anc)

    return pos

def set_attocube_output(axis, state, anc=None):
    '''
    ground or unground attocube
    '''

    anc_passed = True
    if anc==None:
        anc = initialize_attocube()
        anc_passed = False
    
    anc.setAxisOutput(axis, state, 0)

    if anc_passed == False:
        print(anc.getPosition(axis_dict[axis])*1e6)
        close_attocube(anc)

def initialize_attocube():
    '''
    attocube initialization function, which checks if an initialization already exists at attocube_handle before initalizizing another instance.
    '''

    try:
        if os.path.exists(ATTOCUBE_HANDLE_FNAME):
            with open(ATTOCUBE_HANDLE_FNAME, 'rb') as f:
                anc = pickle.load(f)
                read_attocube('x',anc) # test the connection
        else:
            anc = Positioner()
            with open(ATTOCUBE_HANDLE_FNAME, 'wb') as f:
                pickle.dump(anc, f)
    except Exception as e:
        #print(f'initialization error: {e}')
        os.remove(ATTOCUBE_HANDLE_FNAME)
        anc = Positioner()
        with open(ATTOCUBE_HANDLE_FNAME, 'wb') as f:
            pickle.dump(anc, f)

    return anc

def close_attocube(anc):
        if os.path.exists(ATTOCUBE_HANDLE_FNAME):
            os.remove(ATTOCUBE_HANDLE_FNAME)
            anc.close()

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

def move_x(position,  anc=None, tolerance=1, go_back=0):
    '''
    wrapper to move attocube x positioner.
    '''
    move_attocube('x', position, anc, tolerance, go_back)

def move_y(position,  anc=None, tolerance=1, go_back=0):
    '''
    wrapper to move attocube y positioner.
    '''
    move_attocube('y', position, anc, tolerance, go_back)

def move_z(position,  anc=None, tolerance=1, go_back=0):
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
