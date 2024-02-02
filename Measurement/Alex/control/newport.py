import instruments.newport as newport
from pymeasure.instruments.newport import ESP300
from OrensteinLab_git.configuration import config_dict
import time

esp_port = config_dict['ESP COM Port']
esp_model = config_dict['ESP Model']

######################
### Core Functions ###
######################

def move_axis(axis_index, pos, axis=None):
    if esp_model=='301':
        move_esp301(axis_index, pos, axis)
    elif esp_model=='300':
        move_esp300(axis_index, pos, axis)
    else:
        raise ValueError(f'ESP Model {esp_model} not supported.')

def read_axis(axis_index, axis=None, print_flag=True):
    if esp_model=='301':
        pos = read_esp301(axis_index, axis, print_flag)
    elif esp_model=='300':
        pos = read_esp300(axis_index, axis, print_flag)
    else:
        raise ValueError(f'ESP Model {esp_model} not supported.')
    return pos

def corotate_axes(axis_1_index, axis_2_index, angle_1, angle_2, axis_1=None, axis_2=None):
    if esp_model=='301':
        corotate_axes301(axis_1_index, axis_2_index, angle_1, angle_2, axis_1, axis_2)
    elif esp_model=='300':
        corotate_axes300(axis_1_index, axis_2_index, angle_1, angle_2, axis_1, axis_2)
    else:
        raise ValueError(f'ESP Model {esp_model} not supported.')

def initialize_axis(axis_index):
    if esp_model=='301':
        axis = initialize_esp301(axis_index)
    elif esp_model=='300':
        axis = initialize_esp300(axis_index)
    else:
        raise ValueError(f'ESP Model {esp_model} not supported.')
    return axis

def close_axis(axis):
    return 0

def move_esp301(axis_index, pos, axis=None):
    '''
    rotate waveplate. I can now add any checks that have to happen typically.
    '''
    # initialize axis
    if axis==None:
        axis = initialize_esp301(axis_index)

    while True:
        try:
            axis.move(pos, absolute=True)
            check_axis_stability301(axis, axis_index)
            break
        except:
            axis = initialize_esp301(axis_index)
            print('failed to move axis, trying again')

def read_esp301(axis_index, axis=None, print_flag=True):
    '''
    read angle on an axis
    '''
    # initialize axis
    obj_passed = True
    if axis==None:
        axis = initialize_esp301(axis_index)
        obj_passed = False

    while True:
        try:
            pos = float(axis.position)
            break
        except:
            axis = initialize_esp301(axis_index)
            print('failed to read axis, trying again.')
            #pass
        time.sleep(0.1)

    if print_flag==True and obj_passed==False:
        print(pos)
    return pos

def check_axis_stability301(axis, axis_index):
    '''
    helper function for rotate_axis.
    '''
    while True:
        try:
            if axis.is_motion_done==True:
                break
        except:
            axis = initialize_esp301(axis_index)
            print('failed to check axis stability, trying agian.')
            #pass
        time.sleep(0.1)

def move_esp300(axis_index, pos, axis=None):

    if axis==None:
        axis = initialize_esp300(axis_index)

    while True:
        try:
            axis.position = pos
            check_axis_stability300(axis, axis_index)
            break
        except:
            axis = initialize_esp300(axis_index)
            print('failed to move axis, trying again.')
        time.sleep(0.1)

def read_esp300(axis_index, axis=None, print_flag=True):
    obj_passed = True
    if axis==None:
        axis = initialize_esp300(axis_index)
        obj_passed = False

    while True:
        try:
            pos = float(axis.position)
            break
        except:
            axis = initialize_esp300(axis_index)
            print('failed to read axis, trying again.')
            #pass
        time.sleep(0.1)

    if print_flag==True and obj_passed==False:
        print(pos)
    return pos

def check_axis_stability300(axis, axis_index):
    '''
    helper function for rotate_axis.
    '''
    while True:
        time.sleep(0.1)
        try:
            if axis.motion_done==True:
                break
        except:
            axis = initialize_esp300(axis_index)
            print('failed to check axis stability, trying agian.')
            #pass

def corotate_axes301(axis_1_index, axis_2_index, angle_1, angle_2, axis_1=None, axis_2=None):

    if axis_1==None:
        axis_1 = initialize_esp301(axis_1_index)
    if axis_2==None:
        axis_2 = initialize_esp301(axis_2_index)

    while True:
        try:
            axis_1.move(angle_1, absolute=True)
            axis_2.move(angle_2, absolute=True)
            break
        except:
            print('failed to check axis stability, trying again')
            time.sleep(0.1)
    while True:
        time.sleep(0.1)
        try:
            if axis_1.is_motion_done==True and axis_2.is_motion_done==True:
                break
        except:
            axis_1 = initialize_esp301(axis_1_index)
            axis_2 = initialize_esp301(axis_2_index)
            print('failed to check axis stability, trying agian.')
            #pass

def corotate_axes300(axis_1_index, axis_2_index, angle_1, angle_2, axis_1=None, axis_2=None):

    if axis_1==None:
        axis_1 = initialize_esp300(axis_1_index)
    if axis_2==None:
        axis_2 = initialize_esp300(axis_2_index)

        while True:
            try:
                axis_1.position = angle_1
                axis_2.position = angle_2
                break
            except:
                print('failed to check axis stability, trying again')
                time.sleep(0.1)
    while True:
        time.sleep(0.1)
        try:
            if axis_1.motion_done==True and axis_2.motion_done==True:
                break
        except:
            axis_1 = initialize_esp300(axis_1_index)
            axis_2 = initialize_esp300(axis_2_index)
            print('failed to check axis stability, trying agian.')
            #pass

def initialize_esp301(axis_index):

    controller = newport.NewportESP301.open_serial(port=esp_port, baud=921600)
    while True:
        try:
            axis = newport.NewportESP301Axis(controller, axis_index-1)
            axis.enable()
            break
        except:
            time.sleep(0.1)
    return axis

def initialize_esp300(axis_index):

    obj = ESP300(esp_port)
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


#########################
### Wrapper Functions ###
#########################

def move_axis_1(pos, axis=None):
    move_axis(1, pos, axis)

def move_axis_2(pos, axis=None):
    move_axis(2, pos, axis)

def move_axis_3(pos, axis=None):
    move_axis(3, pos, axis)

def corotate_axes12(angle, axes=None, bal_angle=0):
    if axes==None:
        axis_1, axis_2 = initialize_corotate_axes12()
    else:
        axis_1, axis_2 = axes
    corotate_axes(1, 2, angle, angle+bal_angle, axis_1=axis_1, axis_2=axis_2)

def read_axis_1(axis=None, print_flag=True):
    return read_axis(1, axis, print_flag)

def read_axis_2(axis=None, print_flag=True):
    return read_axis(2, axis, print_flag)

def read_axis_3(axis=None, print_flag=True):
    return read_axis(3, axis, print_flag)

def read_corotate_axes12(axes=None, print_flag=True):
    if axes==None:
        axis_1, axis_2 = initialize_corotate_axes12()
    else:
        axis_1, axis_2 = axes
    return read_axis(1, axis_1, print_flag=print_flag)

def initialize_axis_1():
    return initialize_axis(1)

def initialize_axis_2():
    return initialize_axis(2)

def initialize_axis_3():
    return initialize_axis(3)

def initialize_corotate_axes12():
    axis_1 = initialize_axis_1()
    axis_2 = initialize_axis_2()
    return axis_1, axis_2

def close_axis_1(axis):
    return 0

def close_axis_2(axis):
    return 0

def close_axis_3(axis):
    return 0

def close_corotate_axes12(axis):
    return 0
