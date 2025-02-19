import numpy as np
from pyanc350.v4 import Positioner
from OrensteinLab_git.configuration import CONFIG_DICT
import pickle
import time
import os

axis_dict = {'x':0, 'y':1, 'z':2}
ATTOCUBE_HANDLE_FNAME = CONFIG_DICT['Attocube Handle']

######################
### Core Functions ###
######################

def move_attocube(axis, position, anc=None, tolerance=1, go_back=0, ground=False, check_stability=True):
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

    # check stability
    if check_stability:
        while (error >= tol):
            anc.moveAbsolute(axis_dict[axis], int(pos/1e6))
            time.sleep(0.1)
            error = np.abs(pos-anc.getPosition(axis_dict[axis])*1e6)

    # print and close only if another process hasn't passed anc object
    if anc_passed == False:
        print(anc.getPosition(axis_dict[axis])*1e6)
        close_attocube(anc)
    else:
        return anc

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
    else:
        return pos, anc

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

    anc = initialize_attocube()
    move_x(x, anc)
    move_y(y, anc)
    if z==0:
        pass
    else:
        move_z(z, anc)
    close_attocube(anc)
    
def move_x(position,  anc=None, tolerance=1, go_back=0, check_stability=True):
    '''
    wrapper to move attocube x positioner.
    '''
    return move_attocube('x', position, anc, tolerance, go_back, check_stability)

def move_y(position,  anc=None, tolerance=1, go_back=0, check_stability=True):
    '''
    wrapper to move attocube y positioner.
    '''
    return move_attocube('y', position, anc, tolerance, go_back, check_stability)

def move_z(position,  anc=None, tolerance=1, go_back=0, check_stability=True):
    '''
    wrapper to move attocube z positioner.
    '''
    return move_attocube('z', position, anc, tolerance, go_back, check_stability)

def read_x(anc=None, print_flag=True):
    return read_attocube('x', anc, print_flag)

def read_y(anc=None, print_flag=True):
    return read_attocube('y', anc, print_flag)

def read_z(anc=None, print_flag=True):
    return read_attocube('z', anc, print_flag)
