import instruments.newport as newport
from pymeasure.instruments.newport import ESP300
from OrensteinLab_git.configuration import CONFIG_DICT
import time

esp_port = CONFIG_DICT['ESP COM Port']
esp_model = CONFIG_DICT['ESP Model']

######################
### Core Functions ###
######################

def move_axis(axis_index, pos, axis=None, wait_time=0, check_stability=True):

    if esp_model=='301':
        axis =  move_esp301(axis_index, pos, axis, check_stability)
    elif esp_model=='300':
        axis = move_esp300(axis_index, pos, axis, check_stability)
    else:
        raise ValueError(f'ESP Model {esp_model} not supported.')
    time.sleep(wait_time)

    return axis

def read_axis(axis_index, axis=None):

    if esp_model=='301':
        pos, axis = read_esp301(axis_index, axis)
    elif esp_model=='300':
        pos, axis = read_esp300(axis_index, axis)
    else:
        raise ValueError(f'ESP Model {esp_model} not supported.')

    return pos, axis

def corotate_axes(axis_1_index, axis_2_index, angle_1, angle_2, axis_1=None, axis_2=None, check_stability=True):
    if esp_model=='301':
        axis_1, axis_2 = corotate_axes301(axis_1_index, axis_2_index, angle_1, angle_2, axis_1, axis_2, check_stability)
    elif esp_model=='300':
        axis_1, axis_2 = corotate_axes300(axis_1_index, axis_2_index, angle_1, angle_2, axis_1, axis_2, check_stability)
    else:
        raise ValueError(f'ESP Model {esp_model} not supported.')

    return axis_1, axis_2

def initialize_axis(axis_index):
    if esp_model=='301':
        axis = initialize_esp301(axis_index)
    elif esp_model=='300':
        axis = initialize_esp300(axis_index)
    else:
        raise ValueError(f'ESP Model {esp_model} not supported.')
    return axis

def close_axis(axis):
    if esp_model=='301':
        close_esp301(axis)
    elif esp_model=='300':
        close_esp300(axis)
    else:
        raise ValueError(f'ESP Model {esp_model} not supported.')

def move_esp301(axis_index, pos, axis=None, check_stability=True):
    '''
    rotate waveplate. I can now add any checks that have to happen typically.
    '''
    # initialize axis
    axis_passed = True
    if axis==None:
        axis = initialize_esp301(axis_index)
        axis_passed=False

    while True:
        try:
            axis.move(pos, absolute=True)
            if check_stability:
                axis = check_axis_stability301(axis, axis_index)
            break
        except:
            print(f'failed to move axis {axis_index}, trying again')
            close_esp301(axis)
            axis = initialize_esp301(axis_index)
            print(f'reinitialized axis {axis_index}')
            if check_stability:
                axis = check_axis_stability301(axis, axis_index)
            break
        time.sleep(0.1)

    if axis_passed==False:
        close_esp301(axis)
        return None
    else:
        return axis

def read_esp301(axis_index, axis=None):
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
            print(f'failed to read axis {axis_index}, trying again')
            close_esp301(axis)
            axis = initialize_esp301(axis_index)
            #print(f'reinitialized axis {axis_index} at {axis}')
            #pass
        time.sleep(0.1)

    if obj_passed==False:
        close_esp300(axis)
        return pos, None
    else:
        return pos, axis

def check_axis_stability301(axis, axis_index):
    '''
    helper function for rotate_axis.
    '''
    while True:
        try:
            if axis.is_motion_done==True:
                break
        except:
            print(f'failed to check stability axis {axis_index}, trying again')
            close_esp301(axis)
            axis = initialize_esp301(axis_index)
            #print(f'reinitialized axis {axis_index}')
            #pass
        time.sleep(0.1)
    return axis

def move_esp300(axis_index, pos, axis=None, check_stability=True):

    axis_passed = True
    if axis==None:
        axis = initialize_esp300(axis_index)
        axis_passed=False

    while True:
        try:
            axis.position = pos
            if check_stability:
                axis = check_axis_stability300(axis, axis_index)
            break
        except:
            print(f'failed to move axis {axis_index}, trying again')
            close_esp300(axis)
            axis = initialize_esp300(axis_index)
            if check_stability:
                axis = check_axis_stability300(axis, axis_index)
            #print(f'reinitialized axis {axis_index}')
        time.sleep(0.1)

    if axis_passed==False:
        close_esp300(axis)
        return None
    else:
        return axis

def read_esp300(axis_index, axis=None):
    obj_passed = True
    if axis==None:
        axis = initialize_esp300(axis_index)
        obj_passed = False

    while True:
        try:
            pos = float(axis.position)
            break
        except:
            print(f'failed to read axis {axis_index}, trying again')
            close_esp300(axis)
            axis = initialize_esp300(axis_index)
            #print(f'reinitialized axis {axis_index}')
            #pass
        time.sleep(0.1)

    if obj_passed==False:
        close_esp300(axis)
        return pos, None
    else:
        return pos, axis

def check_axis_stability300(axis, axis_index):
    '''
    helper function for rotate_axis.
    '''
    while True:
        try:
            if axis.motion_done==True:
                break
        except:
            print(f'failed to check stability axis {axis_index}, trying again')
            close_esp300(axis)
            axis = initialize_esp300(axis_index)
            #print(f'reinitialized axis {axis_index}')
        time.sleep(0.1)
            #pass
    return axis

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
            print('failed to corotate axes, trying agian.')
            close_esp301(axis_1)
            close_esp301(axis_2)
            axis_1 = initialize_esp301(axis_1_index)
            axis_2 = initialize_esp301(axis_2_index)
            #print('reinitialized axes')
        time.sleep(0.1)
    while True:
        try:
            if axis_1.is_motion_done==True and axis_2.is_motion_done==True:
                break
        except:
            print('failed to check axis corotate axes stability, trying agian.')
            close_esp301(axis_1)
            close_esp301(axis_2)
            axis_1 = initialize_esp301(axis_1_index)
            axis_2 = initialize_esp301(axis_2_index)
            print('reinitialized axes')
        time.sleep(0.1)
            #pass

    return axis_1, axis_2

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
            print('failed to corotate axes, trying agian.')
            close_esp300(axis_1)
            close_esp300(axis_2)
            axis_1 = initialize_esp300(axis_1_index)
            axis_2 = initialize_esp300(axis_2_index)
            #print('reinitialized axes')
        time.sleep(0.1)

    while True:
        try:
            if axis_1.motion_done==True and axis_2.motion_done==True:
                break
        except:
            print('failed to check axis corotate axes stability, trying agian.')
            close_esp300(axis_1)
            close_esp300(axis_2)
            axis_1 = initialize_esp300(axis_1_index)
            axis_2 = initialize_esp300(axis_2_index)
            print('reinitialized axes')
        time.sleep(0.1)
            #pass

    return axis_1, axis_2

def initialize_esp301(axis_index):

    created_controller = False
    while True:
        try:
            #print('connection to controller')
            controller = newport.NewportESP301.open_serial(port=esp_port, baud=921600)
            #print('creating axis')
            created_controller = True
            axis = newport.NewportESP301Axis(controller, axis_index-1)
            #print('enabling axis')
            axis.enable()
            #print('breaking')
            break
        except:
            if created_controller==True:
                controller._file._conn.close()
            print(f'failed to initialize axis {axis_index}, trying again')
            time.sleep(0.1)
    return axis

def initialize_esp300(axis_index):

    while True:
        try:
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
            break
        except:
            print(f'failed to initialize axis {axis_index}, trying again')
            time.sleep(0.1)
    return axis

def close_esp301(axis):
    axis._controller._file._conn.close()

def close_esp300(axis):
    return 0

#########################
### Wrapper Functions ###
#########################

def move_axis_1(pos, axis=None, wait_time=0, check_stability=True):
    return move_axis(1, pos, axis, wait_time, check_stability)

def move_axis_2(pos, axis=None, wait_time=0, check_stability=True):
    return move_axis(2, pos, axis, wait_time, check_stability)

def move_axis_3(pos, axis=None, wait_time=0, check_stability=True):
    return move_axis(3, pos, axis, wait_time, check_stability)

def corotate_axes12(angle, axes=None, bal_angle=0, check_stability=True):
    if axes==None:
        axis_1, axis_2 = initialize_corotate_axes12()
    else:
        axis_1, axis_2 = axes
    return corotate_axes(1, 2, angle, angle+bal_angle, axis_1=axis_1, axis_2=axis_2, check_stability=check_stability)

def read_axis_1(axis=None):
    return read_axis(1, axis)

def read_axis_2(axis=None):
    return read_axis(2, axis)

def read_axis_3(axis=None):
    return read_axis(3, axis)

def read_corotate_axes12(axes=None):
    if axes==None:
        axis_1, axis_2 = initialize_corotate_axes12()
    else:
        axis_1, axis_2 = axes
    pos, axis_1 = read_axis(1, axis_1)
    return pos, [axis_1, axis_2]

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
    close_axis(axis)

def close_axis_2(axis):
    close_axis(axis)

def close_axis_3(axis):
    close_axis(axis)

def close_corotate_axes12(axis):
    ax1, ax2 = axis
    close_axis(ax1)
    close_axis(ax2)
