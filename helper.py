'''
helper methods
'''

######################
### Helper Methods ###
######################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import time
import threading
from tqdm.auto import tqdm
import scipy.optimize as opt
import scipy.interpolate as interp
import inspect
import pickle
from OrensteinLab_git.devices.concurrency_classes import LockedVar, StoppableThread, LockedDict
from OrensteinLab_git.motors_and_instruments import MOTOR_DICT, INSTRUMENT_DICT, META_MOTORS, ACTIVE_MOTORS, ACTIVE_INSTRUMENTS

###
### Instrument Methods
###

def capture_instrument_information(inst_dict):
    '''
    creates ikwargs_dict from inst_dict and ACTIVE_INSTRUMENTS

    args:
        - inst_dict:        dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}}

    returns:
        - ikwargs_dict:     dictionary of kwarg_dicts for the read function of each instrument, ex) {'lockin':{'time_constant':0.3}}
    '''
    
    instrument_names = list(inst_dict.keys())
    ikwargs_dict={}
    for i in ACTIVE_INSTRUMENTS:
        if i in instrument_names:
            ikwargs_dict[i]=inst_dict[i]
        else:
            ikwargs_dict[i]={}

    return ikwargs_dict

def initialize_instruments(iobj_dict):
    '''
    given iobj_dict with instrument objects that have already been initialized, returns iobj_dict with all instruments in ACTIVE_INSTRUMENTS initialized

    args:
        - iobj_dict:        dictionary of instrument handle objects

    returns:
        - iobj_dict:        modified dictionary of instrument handle objects including all instruments in ACTIVE_INSTRUMENTS
    '''

    init_instruments = list(iobj_dict.keys())
    iobj_dict_new={}
    for i in ACTIVE_INSTRUMENTS:
        if i in init_instruments:
            iobj_dict_new[i] = iobj_dict[i]
        else:
            init_func = INSTRUMENT_DICT[i]['init']
            iobj_dict_new[i] = init_func()
    return iobj_dict_new

def get_instruments_header(iobj_dict, ikwargs_dict):
    '''
    given an instrument object dictionary, returns list of strings which are used to generate a file header of all the measured values of the instrument

    args:
        - iobj_dict:        dictionary of instrument handle objects
        - ikwargs_dict:     dictionary of kwarg_dicts for the read function of each instrument, ex) {'lockin':{'time_constant':0.3}}

    returns:
        - instrument_header:    list of variable names to be read by instruments
    '''
    data_dict = read_instruments(iobj_dict, ikwargs_dict)
    instruments_header = list(data_dict.keys())
    return instruments_header

def read_instruments(iobj_dict, ikwargs_dict):
    '''
    given iobj_dict, reads all instruments with kwargs given in iobj_dict and returns read values as a dictionary.
    NOTE: it is the user's responsibility to ensure that reading instruments never yields repeated names for a measured instrument

    args:
        - iobj_dict:        dictionary of instrument handle objects
        - ikwargs_dict:     dictionary of kwarg_dicts for the read function of each instrument, ex) {'lockin':{'time_constant':0.3}}

    returns:
        - data_dict:        dictionary of data, where keys are variable names and values are variables measured values
    '''
    instruments = list(iobj_dict.keys())
    data_dict={}
    for i in instruments:
        iobj = iobj_dict[i]
        read_func = INSTRUMENT_DICT[i]['read']
        instrument_data = read_func(iobj, **ikwargs_dict[i])
        for d in list(instrument_data.keys()):
            data_dict[d] = instrument_data[d]
    return data_dict

def close_instruments(iobj_dict):
    '''
    given a dictionary of instruments objects, closes all instruments=list(iobj_dict.keys())

    args:
        - iobj_dict:        dictionary of instrument handle objects
    '''
    instruments = list(iobj_dict.keys())
    for i in instruments:
        obj = iobj_dict[i]
        close_func = INSTRUMENT_DICT[i]['close']
        close_func(obj)

###
### Motor Methods
###

def get_motor_range(start, end, step_size):
    '''
    returns a numpy array of positions corresponding to input. The convention is that step_size is always a positive number and that the motor moves from start to end.

    args:
        - start:        starting position
        - end:          ending position
        - step_size     step size
    
    return:
        - range:        numpy array of positions corresponding to arguments
    '''
    step_size=abs(step_size)
    if start < end:
        dir = 1
    else:
        dir = -1
    list1 = np.arange(start, end, dir*step_size)
    list2 = np.array([end])
    range = np.concatenate((list1, list2))
    return range

def capture_motor_information(map_dict):
    '''
    unpack map_dict used to define a motor scan

    args:
        - map_dict:         dictionary of key:value pairs where key is name of a motor in MOTOR_DICT and value is a tuple (start, stop step, kwargs) or (positions, kwargs), and where kwargs is another dictionary where keys are kwarg names for the motor move function and value the corresponding value to set during move steps. motors are scanned from left to right.

            ex: {'temp':(100,200,10,{'wait_time':30})} or {'temp':(np.arange(100,200,10),{'wait_time':30})}

    return:
        - motors:          list of motors in motor scan
        - mranges:         list of positions for each motor
        - mkwargs_dict:    dictionary of kwarg_dicts for the move function of each motor, ex) {'temp':{'wait_time':30}}

    '''
    

    motors = list(map_dict.keys())

    # check validity and setup motor ranges and kwargs
    mranges = []
    mkwargs_dict = {}
    valid_motors = list(MOTOR_DICT.keys())
    for m in motors:
        if m not in valid_motors:
            raise ValueError(f'Invalid motor name. Please select motors from the list {valid_motors}.')
        if len(map_dict[m])==2:
            range = map_dict[m][0]
            kwargs = map_dict[m][1]
        elif len(map_dict[m])==4:
            start = map_dict[m][0]
            end = map_dict[m][1]
            step_size = map_dict[m][2]
            kwargs = map_dict[m][3]
            range = get_motor_range(start, end, step_size)
        mranges.append(range)
        mkwargs_dict[m] = kwargs

    return motors, mranges, mkwargs_dict

def capture_sequence_information(sequence_list):
    '''
    unpack sequence_list where sequence_list=[('step1_motor', step1_target, step1_kwarg_dict, step1_tol),(step2_motor,...),...]

    args:
        - sequence_list:       list of tuples defining a motor sequence. Each tuple must be of the form ('motor', target, kwargs_dict, tolerance), where tolerances are on the target position

    return:
        - sequence_motors:     list of motor that is moved in each step in sequence
        - mtargets:            list of target positions for each step in sequence
        - mkwrags:             list of kwarg dictionaries for each step in sequence. ex) [{'wait_time':30}, {'wait_time':30},...]
        - motors:              list of unique motors moved during sequence
        - tols:                list of tolerances for each step in sequence

    '''
    sequence_motors, mtargets,mkwargs,tols = [[],[],[],[]]
    for move in sequence_list:
        m = move[0]
        valid_motors = list(MOTOR_DICT.keys())
        if m not in valid_motors:
            raise ValueError(f'Invalid motor name. Please select motors from the list {valid_motors}.')
        sequence_motors.append(m)
        mtargets.append(move[1])
        mkwargs.append(move[2])
        tols.append(move[3])
    motors = list(set(sequence_motors))
    return sequence_motors, mtargets, mkwargs, motors, tols

def initialize_motors(motors, mobj_dict):
    '''
    given motor list and mobj_dict of some or all initialized motors, returns mobj_dict with all motors initialized

    args:
        - motors:           list of motors
        - mobj_dict:        dictionary of motor handle object; can include any number of motors

    returns:
        - mobj_dict_new:    dictionary of motor handle objects, including initialized motors in variable motors which was not in mobj_dict to begin with
    '''
    init_motors = list(mobj_dict.keys())
    mobj_dict_new = {}
    for m in motors:
        if m in init_motors:
            mobj_dict_new[m] = mobj_dict[m]
        else:
            init_func = MOTOR_DICT[m]['init']
            mobj_dict_new[m] = init_func()
    return mobj_dict_new

def initialize_motors_threaded(motors):
    '''
    UNDER CONSTRUCTION
    '''
    mobj_dict_locked = LockedDict({})
    initialization_threads = []
    for m in motors:
        init_func = MOTOR_DICT[m]['init']
        init_thread = threading.Thread(target=add_initialized_motor_to_thread_dict, args=(m, init_func, mobj_dict_locked))
        initialization_threads.append(init_thread)
    for ii, init_thread in enumerate(initialization_threads):
        print(f'starting initialization thread for motor {motors[ii]}')
        init_thread.start()
    for ii, init_thread in enumerate(initialization_threads):
        init_thread.join()
        print(f'joined initialization thread for motor {motors[ii]}')
    # sort according to motors
    mobj_dict_unsorted = mobj_dict_locked.locked_get_dict()
    mobj_dict = dict(sorted(mobj_dict_unsorted.items(), key=lambda x:motors.index(x[0])))
    return mobj_dict

def add_initialized_motor_to_thread_dict(motor, init_func, threaded_dict):
    '''
    UNDER CONSTRUCTION
    '''
    while True:
        try:
            print(f'starting initialization of motor {motor}')
            obj = init_func()
            threaded_dict.locked_update(motor, obj)
            print(f'finished initializing motor {motor}')
            break
        except:
            time.sleep(0.1)

def close_motors(mobj_dict):
    '''
    given a dictionary of motor objects, closes all motors=list(mobj_dict.keys())

    args:
        - mobj_dict:        dictionary of motor handle objects
    '''
    motors = list(mobj_dict.keys())
    for m in motors:
        obj = mobj_dict[m]
        close_func = MOTOR_DICT[m]['close']
        close_func(obj)

def move_motors_to_start(motors, mobj_dict, mkwargs_dict, positions, print_flag=False):
    '''
    moves a set of motors to initial positions. A "move back" move is invoked for handling hysteretic motors.

    args:
        - motors:           list of motors
        - mobj_dict:        dictionary of motor handle objects; can include any number of motors
        - mkwrags_dict:     dictionary of kwargs for each motor's move function
        - positions:        list of target positions
        - print_flag:       if true, prints position to which motor has moved
    '''
    for ii, m in enumerate(motors):
        move_back = MOTOR_DICT[m]['move_back']
        p = positions[0][ii]
        move_func = MOTOR_DICT[m]['move']
        read_func = MOTOR_DICT[m]['read']
        obj = mobj_dict[m]
        kwargs = mkwargs_dict[m]
        move_func(p - move_back, obj, **kwargs)
        move_func(p, obj, **kwargs)
        p_true = read_func(obj)
        if print_flag==True:
            print(f'Moved motor {m} to {p_true}.')

def move_motors(motors, mobj_dict, mkwargs_dict, current_pos, start_pos, new_pos, print_flag=False):
    '''
    moves a set of motors. if motors are hysteretic, a "move back" move is invoked if the motor is directed to move to a defined initial position

    args:
        - motors:           list of motors
        - mobj_dict:        dictionary of motor handle objects; can include any number of motors
        - mkwrags_dict:     dictionary of kwargs for each motor's move function
        - current_pos:      list of motor positions before move
        - start_pos:        list of motor starting positions (for invoking "move back" motions for hysteretic motors)
        - new_pos:          list of motor target positions
        - print_flag:       if true, prints position to which motor has moved
    '''
    for ii, m in enumerate(motors):
        p_old = current_pos[ii]
        p_new = new_pos[ii]
        p_start = start_pos[ii]
        if p_new!=p_old:
            move_func = MOTOR_DICT[m]['move']
            read_func = MOTOR_DICT[m]['read']
            obj = mobj_dict[m]
            kwargs = mkwargs_dict[m]
            if p_new==p_start:
                move_back = MOTOR_DICT[m]['move_back']
                move_func(p_new-move_back, obj, **kwargs)
                move_func(p_new, obj, **kwargs)
            else:
                move_func(p_new, obj, **kwargs)
            if print_flag==True:
                p_true = read_func(obj)
                print(f'Moved motor {m} to {p_true}.')

def read_motors(motors, mobj_dict):
    '''
    given a list of motors and mobj_dict of motor objects, returns a dictionary of values for each motor in motors

    args:
        - motors:       list of motors
        - mobj_dict:    dictionary of motor handle objects; can include any number of motors
    
    return:
        - pos_dict:     dictionary of motor position values for motors in motors variable.
    '''
    pos_dict = {}
    for m in motors:
        mobj = mobj_dict[m]
        read_func = MOTOR_DICT[m]['read']
        pos_dict[m] = read_func(mobj)
    return pos_dict

def gen_positions_recurse(range_list, n, pos_list=[], current_pos=None):
    '''    given an empty pos_list, and a range_list, recursively generates a list of positions that span the spacce in range_list. Note that positions are written from first entry in range_list to last.

    args:
        - range_list:       a list of np arrays, where each array is a range of interest.
        - n:                max index of arrays in range_list, needed for recursion
        - post_list:        should be an empty list which the function will append to
        - current_pos:      n+1 dim array that carries around the positions to append for each recursive iteration.

    returns:
        - post_list
    '''
    if n==len(range_list)-1:
        current_pos = np.asarray([i[0] for i in range_list])#np.asarray(range_list)[:,0]
        pos_list = []
    if n>=0:
        for i in range_list[n]:
            current_pos[n] = i
            pos_list = gen_positions_recurse(range_list, n-1, pos_list, current_pos)
    else:
        pos_list.append(np.copy(current_pos))

    return pos_list

###
### File Writing Methods
###

def generate_filename(filename_head, filename, current_func):
    '''
    Generates filename_head and filename. If we supply it with no filename head or filename it autogenerates both. If we only supply a filaname it inferes the head to be current working directory.

    args:
        - filename_head:    path to which file is generated
        - filename:         name of file
        - current_func:     name of current function. used in event that filename=None to autogenerate a filename

    return:
        - filename_head:    new path to which file is generated - autogenerated files get placed in filename_head/autogenerated
        - filename:         new name of file    
    '''
    if filename_head==None:
        if filename==None:
            filename_head = os.path.join(os.getcwd(), 'autogenerated')
            if os.path.isdir(filename_head)==False:
                os.mkdir(filename_head)
            filename = current_func
        else:
            filename_head = os.getcwd()
    else:
        if filename==None:
            filename = current_func
    return filename_head, filename

def get_unique_filename(filename_head, filename):
    '''
    give filename_head and filename, return a unique filename with numbers appended to end if necessary

    args:
        - filename_head:    path to which file is generated
        - filename:         name of file

    return:
        - filename_head:    new path to which file is generated
        - filename:         new name of file   
    '''
    #fname = f'{filename_head}\\{filename}.dat'
    fname = os.path.join(filename_head, f'{filename}.dat')
    if not os.path.exists(fname):
        return fname

    path, name = os.path.split(fname)
    name, ext = os.path.splitext(name)

    make_fn = lambda i: os.path.join(path, '%s_%s%s' % (name, i, ext))

    for i in range(1, 1000):
        uni_fn = make_fn(i)
        if not os.path.exists(uni_fn):
            return uni_fn

def write_file_header(fname, header, metadata=None):
    '''
    Helper function that takes in metadata (a dictionary) and header (a list of strings) and writes data in the form:

    [Metadata]
    key: value
    ...
    [Data]
    header

    header strings are written into file tab separated

    args:
        - fname:        full path of file to write
        - header:       list of string representing data variables
        - metadata:     dictionary containing metadata entries
    '''
    with open(fname, 'w') as file:
        if metadata is not None:
            file.write(f'[Metadata]\n')
            for key in list(metadata.keys()):
                file.write(f'{key}:\t{metadata[key]}\n')
            file.write(f'[Data]\n')
        for h in header[:-1]:
            file.write(f'{h}\t')
        file.write(f'{header[-1]}\n')

def append_data_to_file(fname, values):
    '''
    helper function to add a line to a file given a list of values

    values are written into file tab separated

    args:
        - fname:    full path to file
        - values:   list of values to write into file
    '''
    with open(fname, 'a') as file:
        for val in values[:-1]:
            file.write(f'{format(val, ".15f")}\t')
        file.write(f'{format(values[-1], ".15f")}\n')

def save_data_to_file(fname, data, header, metadata=None):
    '''
    utility function for saving data to a file, with optional metadata

    args:
        - fname(string):           full path to datafile
        - data(array):             (n,m) array containing data
        - header(array):           (m) array of strings labeling each column
        - metadata(dict):          a dictionary of parameters to store at the top of the datafile under [Metadata]

    returns: None
    '''
    if not(len(header) == len(data[0,:])):
        raise ValueError('number of header items does not match number of data columns.')
    with open(fname, 'w') as f:
        if not(metadata==None):
            f.write('[METADATA]\n')
            for key, value in list(metadata.keys()):
                f.write(f'{key}:\t{value}\n')
            f.write('[DATA]\n')
        for item in header:
            f.write(str(item)+'\t')
        f.write('\n')
        for line in data:
            for item in line:
                f.write(str(item)+'\t')
            f.write('\n')