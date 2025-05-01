'''
Methods for orchestrating measurements and lab utility methods
'''
get_ipython().run_line_magic('matplotlib', 'notebook')
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
import OrensteinLab_git.devices.newport as newport
from OrensteinLab_git.concurrency import LockedVar, StoppableThread, LockedDict
from OrensteinLab_git.motors_and_instruments import MOTOR_DICT, INSTRUMENT_DICT, ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, DEFAULT_VARS
import OrensteinLab_git.helper as helper
import OrensteinLab_git.plotters as plotters
import OrensteinLab_git.concurrency as concurrency
from orenstein_analysis.measurement import loader, process

#################################################
###                                           ###
###   General Methods Applicable In Any Lab   ###
###                                           ###
#################################################

############################
### Measurement Medthods ###
############################

def motor_scan(map_dict, mkwargs_read_dict={}, ikwargs_dict={}, mobj_dict={}, iobj_dict={}, vars=[], metadata={}, filename_head=None, filename=None, savefile=True, print_flag=False, plot=True, close_devices=True):
    '''
    method to record measurements on instruments as a function of motors specified by dictionary map_dict. meant to acquire data one point at a time in a multidimensional motor space.

    args:
        - map_dict:         dictionary of key:value pairs where key is name of a motor in MOTOR_DICT and value is a tuple (start, stop step, mkwargs_move) or (positions, mkwargs_move), and where mkwargs_move is a dictionary where keys are kwarg names for the motor move function, and value the corresponding value to set during move calls. motors are scanned from left to right.

            ex: {'temp':(100,200,10,{'wait_time':30})} or {'temp':(np.arange(100,200,10),{'wait_time':30})}

        - mkwargs_read_dict:    dictionary of mkwargs for motor read functions. if left empty, defaults to no kwargs
        - ikwargs_dict:         dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}}

        - mobj_dict:            dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:            same as mobj_dict but for instruments
        - vars:                 list of variable names to return and, if applicable, plot, must match something in measurd_header (see below). defaults to DEFAULT_VARS in configuration.py
        - metadata:             dictionary with any additional entries to add to metadata
        - filename_head:        path to directory in which to save data. defaults ot current working directory of notebook
        - filename:             name of file
        - savefile:             if true, saves file
        - print_flag:           if true, print motor positions during scan
        - plot:                 if true, plots data in vars
        - close_devices:        if true, closes all devices in mobj_dict and iobj_dict
    '''

    # parse motor scanning information contained in map_dict
    motors, mranges, mkwargs_move_dict = helper.capture_motor_information(map_dict)

    # initialize motors and instruments, and create kwarg dictionaries for each
    mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict = helper.initialize_instruments_and_motors(ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict)

    # create motor and instrument header
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
    motors_header = [MOTOR_DICT[m]['name'] for m in motors] # for tracking idealized motor positionn
    motors_header_measured = [MOTOR_DICT[m]['name']+str(' measured') for m in ACTIVE_MOTORS] # for tracking measured motor positions
    header = motors_header+instruments_header+motors_header_measured

    # take default value for vars, and check that vars is valid
    if vars==[]:
        vars = DEFAULT_VARS
    if not set(vars).issubset(set(header)):
        raise ValueError(f'vars {vars} not in measurd values {header}')

    # setup file for writing - adds metadata at top and writes the data header
    if savefile:
        # setup metadata
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        metamotors_dict = dict([(MOTOR_DICT[m]['name'], measured_motors_dict[m]) for m in list(measured_motors_dict.keys())])
        metadata = {**metadata, **metamotors_dict}
        # setup and write metadta + header to file
        fname = helper.setup_filename(filename, filename_head, 'motor_scan_'+'_'.join([f'{m}' for m in motors]))
        helper.write_file_header(fname, header, metadata)

    # setup plots for displaying
    if plot==True:
        nvars = len(vars)
        if len(map_dict)==1: # setup 1 dimensional plots for single motor scans
            xlabel = MOTOR_DICT[motors[0]]['name']
            xrange = mranges[0]
            fig, axes, plot_handles_dict, vdata1d_dict = plotters.setup_1d_plots(vars, xlabel, xrange)
        elif len(map_dict)==2: # setup 2 dimensional plots for double motor scans
            xrange = mranges[0]
            yrange = mranges[1]
            xlabel = MOTOR_DICT[motors[0]]['name']
            ylabel = MOTOR_DICT[motors[1]]['name']
            xn = len(xrange)
            fig, axes, plot_handles_dict, vdata2d_dict = plotters.setup_2d_plots(vars, xlabel, xrange, ylabel, yrange)
        else:
            print('Cannot plot scans that are greater than 2 dimensions.')

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    # generate motor scan positions recursively.
    # positions: list of positions where each element contains positions of each motor for nth step in scan
    positions = helper.gen_positions_recurse(mranges, len(mranges)-1)
    num_pos = len(positions)

    # move motors to start
    move_dict = {}
    for jj, m in enumerate(motors):
        move_dict[m] = (positions[0][jj], True)
    mobj_dict = helper.move_motors(move_dict, mobj_dict, mkwargs_move_dict, mkwargs_read_dict, print_flag=print_flag)

    # loop over positions, only moving a motor if its target position has changed.
    start_pos = positions[0]
    current_pos = start_pos
    for ii in tqdm(range(num_pos)):
        pos = positions[ii]

        # if run is False, break loop
        if not run.locked_read():
            break

        # move motors for which position has changed, moving first to "go back" position if moving back to start of initial motor position
        move_dict = {}
        for jj, m in enumerate(motors):
            mtarget = pos[jj]
            if current_pos[jj]!=mtarget:
                move_back_flag=False
                if start_pos[jj]==mtarget:
                    move_back_flag=True
                move_dict[m] = (mtarget, move_back_flag)
        mobj_dict = helper.move_motors(move_dict, mobj_dict, mkwargs_move_dict, mkwargs_read_dict, print_flag=print_flag)
        current_pos = pos

        # acquire data and read actual motor positions - should test how long this takes for typical motor setup
        #t0 = time.time()
        instrument_data_dict, iobj_dict = helper.read_instruments(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        #tf = time.time()
        #print(tf-t0)

        # combine dictionaries to keep all new data in a single dictionary
        newdata = {**instrument_data_dict, **measured_motors_dict}

        # convert acquired data from dictionaries to list for file writing
        instrument_data = [instrument_data_dict[h] for h in instruments_header]
        measured_motors_data = [measured_motors_dict[m] for m in ACTIVE_MOTORS]

        # append new data to file
        if savefile:
            helper.append_data_to_file(fname, list(pos)+instrument_data+measured_motors_data)

        # update plots
        if plot==True:
            if len(map_dict)==1:
                xind = ii
                vdata1d_dict = plotters.update_1d_plots(fig, axes, vars, plot_handles_dict, vdata1d_dict, xind, newdata)
            elif len(map_dict)==2:
                xind = int(ii%xn)
                yind = int(ii/xn)
                vdata2d_dict = plotters.update_2d_plots(fig, axes, vars, plot_handles_dict, vdata2d_dict, xind, yind, newdata)

    # wait for thread to finish
    run.locked_update(False)
    user_input_thread.join()

    # close instruments and motors
    if close_devices:
        helper.close_instruments(iobj_dict)
        helper.close_motors(mobj_dict)
    else:
        return iobj_dict, mobj_dict

def motorfunc_scan(map_dict, mkwargs_read_dict={}, ikwargs_dict={}, mobj_dict={}, iobj_dict={}, vars=[], metadata={}, filename_head=None, filename=None, savefile=True, print_flag=False, plot=True, close_devices=True, function=None):
    '''
    method to record measurements on instruments as a function of motors specified by dictionary map_dict. meant to acquire data one point at a time in a multidimensional motor space.

    very similar to motor_scan, except acts 'function' after moving motors at each step. function can be any use defined python function which has arguments (position, mobj_dict, iobj_dict, mkwrags_read_dict, ikwargs_dict, mkwargs_move_dict), where position is list of current motor positions and returns the same arguments.

    args:
        - map_dict:         dictionary of key:value pairs where key is name of a motor in MOTOR_DICT and value is a tuple (start, stop step, mkwargs_move) or (positions, mkwargs_move), and where mkwargs_move is a dictionary where keys are kwarg names for the motor move function, and value the corresponding value to set during move calls. motors are scanned from left to right.

            ex: {'temp':(100,200,10,{'wait_time':30})} or {'temp':(np.arange(100,200,10),{'wait_time':30})}

        - mkwargs_read_dict:    dictionary of mkwargs for motor read functions. if left empty, defaults to no kwargs
        - ikwargs_dict:         dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}

        - mobj_dict:            dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:            same as mobj_dict but for instruments
        - vars:                 list of variable names to return and, if applicable, plot, must match something in measurd_header (see below). defaults to DEFAULT_VARS in configuration.py
        - metadata:             dictionary with any additional entries to add to metadata
        - filename_head:        path to directory in which to save data. defaults ot current working directory of notebook
        - filename:             name of file
        - savefile:             if true, saves file
        - print_flag:           if true, print motor positions during scan
        - plot:                 if true, plots data in vars
        - close_devices:        if true, closes all devices in mobj_dict and iobj_dict
    '''

    # define default function, which does nothing
    if function is None:
        def function(pos, mobj_d, iobj_d):
            return mobj_d, iobj_d

    # parse motor scanning information contained in map_dict
    motors, mranges, mkwargs_move_dict = helper.capture_motor_information(map_dict)

    # initialize motors and instruments, and create kwarg dictionaries for each
    mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict = helper.initialize_instruments_and_motors(ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict)

    # create motor and instrument header
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
    motors_header = [MOTOR_DICT[m]['name'] for m in motors] # for tracking idealized motor positionn
    motors_header_measured = [MOTOR_DICT[m]['name']+str(' measured') for m in ACTIVE_MOTORS] # for tracking measured motor positions
    header = motors_header+instruments_header+motors_header_measured

    # take default value for vars, and check that vars is valid
    if vars==[]:
        vars = DEFAULT_VARS
    if not set(vars).issubset(set(header)):
        raise ValueError(f'vars {vars} not in measurd values {header}')

    # setup file for writing - adds metadata at top and writes the data header
    if savefile:
        # setup metadata
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        metamotors_dict = dict([(MOTOR_DICT[m]['name'], measured_motors_dict[m]) for m in list(measured_motors_dict.keys())])
        metadata = {**metadata, **metamotors_dict}
        # setup and write metadta + header to file
        fname = helper.setup_filename(filename, filename_head, 'motorfunc_scan_'+'_'.join([f'{m}' for m in motors]))
        helper.write_file_header(fname, header, metadata)

    # setup plots for displaying
    if plot==True:
        nvars = len(vars)
        if len(map_dict)==1: # setup 1 dimensional plots for single motor scans
            xlabel = MOTOR_DICT[motors[0]]['name']
            xrange = mranges[0]
            fig, axes, plot_handles_dict, vdata1d_dict = plotters.setup_1d_plots(vars, xlabel, xrange)
        elif len(map_dict)==2: # setup 2 dimensional plots for double motor scans
            xrange = mranges[0]
            yrange = mranges[1]
            xlabel = MOTOR_DICT[motors[0]]['name']
            ylabel = MOTOR_DICT[motors[1]]['name']
            xn = len(xrange)
            fig, axes, plot_handles_dict, vdata2d_dict = plotters.setup_2d_plots(vars, xlabel, xrange, ylabel, yrange)
        else:
            print('Cannot plot scans that are greater than 2 dimensions.')

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    # generate motor scan positions recursively.
    # positions: list of positions where each element contains positions of each motor for nth step in scan
    positions = helper.gen_positions_recurse(mranges, len(mranges)-1)
    num_pos = len(positions)

    # move motors to start
    move_dict = {}
    for jj, m in enumerate(motors):
        move_dict[m] = (positions[0][jj], True)
    mobj_dict = helper.move_motors(move_dict, mobj_dict, mkwargs_move_dict, mkwargs_read_dict, print_flag=print_flag)

    # loop over positions, only moving a motor if its target position has changed.
    start_pos = positions[0]
    current_pos = start_pos
    for ii in tqdm(range(num_pos)):
        pos = positions[ii]

        # if run is False, break loop
        if not run.locked_read():
            break

        # carry out function step before moving and data acquisition
        mobj_dict, iobj_dict, mkwargs_read_dict, ikwargs_dict, mkwargs_move_dict = function(pos, mobj_dict, iobj_dict, mkwargs_read_dict, ikwargs_dict, mkwargs_move_dict)

        # move motors for which position has changed, moving first to "go back" position if moving back to start of initial motor position
        move_dict = {}
        for jj, m in enumerate(motors):
            mtarget = pos[jj]
            if current_pos[jj]!=mtarget:
                move_back_flag=False
                if start_pos[jj]==mtarget:
                    move_back_flag=True
                move_dict[m] = (mtarget, move_back_flag)
        mobj_dict = helper.move_motors(move_dict, mobj_dict, mkwargs_move_dict, mkwargs_read_dict, print_flag=print_flag)
        current_pos = pos

        # acquire data and read actual motor positions - should test how long this takes for typical motor setup
        #t0 = time.time()
        instrument_data_dict, iobj_dict = helper.read_instruments(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        #tf = time.time()
        #print(tf-t0)

        # combine dictionaries to keep all new data in a single dictionary
        newdata = {**instrument_data_dict, **measured_motors_dict}

        # convert acquired data from dictionaries to list for file writing
        instrument_data = [instrument_data_dict[h] for h in instruments_header]
        measured_motors_data = [measured_motors_dict[m] for m in ACTIVE_MOTORS]

        # append new data to file
        if savefile:
            helper.append_data_to_file(fname, list(pos)+instrument_data+measured_motors_data)

        # update plots
        if plot==True:
            if len(map_dict)==1:
                xind = ii
                vdata1d_dict = plotters.update_1d_plots(fig, axes, vars, plot_handles_dict, vdata1d_dict, xind, newdata)
            elif len(map_dict)==2:
                xind = int(ii%xn)
                yind = int(ii/xn)
                vdata2d_dict = plotters.update_2d_plots(fig, axes, vars, plot_handles_dict, vdata2d_dict, xind, yind, newdata)

    # wait for thread to finish
    run.locked_update(False)
    user_input_thread.join()

    # close instruments and motors
    if close_devices:
        helper.close_instruments(iobj_dict)
        helper.close_motors(mobj_dict)
    else:
        return iobj_dict, mobj_dict

def timed_measurement(recording_time, mkwargs_read_dict={}, ikwargs_dict={}, mobj_dict={}, iobj_dict={}, vars=[], metadata={}, xvar='time', filename_head=None, filename=None, savefile=True, plot=True, close_devices=True):
    '''
    makes a measurement for a specified recording time.

    args:
        - recording_time:   time in s to record for
        - mkwargs_read_dict:    dictionary of mkwargs for motor read functions. if left empty, defaults to no kwargs
        - ikwargs_dict:        dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}

        - mobj_dict:            dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:            same as mobj_dict but for instruments
        - vars:                 list of variable names to return and, if applicable, plot, must match something in measurd_header (see below). defaults to DEFAULT_VARS in configuration.py
        - metadata:             dictionary with any additional entries to add to metadata
        - xvar:                 variable to plot on x axis. defaults time
        - filename_head:        path to directory in which to save data. defaults ot current working directory of notebook
        - filename:             name of file
        - savefile:             if true, saves file
        - plot:                 if true, plots data in vars
        - close_devices:        if true, closes all devices in mobj_dict and iobj_dict
    '''

    # initialize motors and instruments, and create kwarg dictionaries for each
    mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict = helper.initialize_instruments_and_motors(ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict)

    # create motor and instrument header
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
    motors_header_measured = [MOTOR_DICT[m]['name']+str(' measured') for m in ACTIVE_MOTORS] # for tracking measured motor positions
    header = ['Time (s)']+instruments_header+motors_header_measured

    # take default value for vars, and check that vars is valid
    if vars==[]:
        vars = DEFAULT_VARS
    if not set(vars).issubset(set(header)):
        raise ValueError(f'vars {vars} not in measurd values {header}')

    # setup file for writing - adds metadata at top and writes the data header
    if savefile:
        # setup metadata
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        metamotors_dict = dict([(MOTOR_DICT[m]['name'], measured_motors_dict[m]) for m in list(measured_motors_dict.keys())])
        metadata = {**metadata, **metamotors_dict}
        # setup and write metadta + header to file
        fname = helper.setup_filename(filename, filename_head, 'timed_measurement')
        helper.write_file_header(fname, header, metadata)

    # setup plots for displaying
    if plot:
        xlabel = 'Time (s)'
        if xvar!='time':
            xlabel = MOTOR_DICT[xvar]['name']
        fig, axes, plot_handles_dict, xrange, vdata1d_dict = plotters.setup_1d_plots_append(vars, xlabel)

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    # loop
    t_delay = 0
    tic = time.perf_counter()
    while (t_delay<recording_time):
        #time.sleep(time_constant)
        toc = time.perf_counter()
        t_delay = toc - tic

        # if run is False, break loop
        if not run.locked_read():
            break

        # acquire data and read actual motor positions - should test how long this takes for typical motor setup
        #t0 = time.time()
        instrument_data_dict, iobj_dict = helper.read_instruments(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        #tf = time.time()
        #print(tf-t0)

        # combine dictionaries to keep all new data in a single dictionary
        newdata = {**instrument_data_dict, **measured_motors_dict}

        # convert acquired data from dictionaries to list for file writing
        instrument_data = [instrument_data_dict[h] for h in instruments_header]
        measured_motors_data = [measured_motors_dict[m] for m in ACTIVE_MOTORS]

        # append new data to file
        if savefile:
            helper.append_data_to_file(fname, list([t_delay])+instrument_data+measured_motors_data)

        # update plots
        if plot==True:
            newx = t_delay
            if xvar!='time':
                newx = newdata[xlabel]
            xrange, vdata1d_dict = plotters.update_1d_plots_append(fig, axes, vars, plot_handles_dict, xrange, vdata1d_dict, newx, newdata)

    # wait for thread to finish
    run.locked_update(False)
    user_input_thread.join()

    # close instruments and motors
    if close_devices:
        helper.close_instruments(iobj_dict)
        helper.close_motors(mobj_dict)
    else:
        return iobj_dict, mobj_dict

def motor_sequence(sequence_list, mkwargs_read_dict={}, ikwargs_dict={}, mobj_dict={}, iobj_dict={}, vars=[], metadata={}, xvar='time', filename_head=None, filename=None, savefile=True, plot=True, close_devices=True):
    '''
    makes a continuos measurement over a sequence of motor steps.

    NOTE: for all motors in sequence, move function must have a check_stability kwargs, which when False induces function to return immediately after setting target position

    args:
        - sequence_list:       list of motor steps to make sequentially, where each entry is a tuple ('motor_name', target, kwargs_move_dict, tolerance)

            ex: [('opticool_field', 1000, {}, 1), ('opticool_field', -1000, {}, 1), ('opticool_field', 0, {}, 1)]

        - mkwargs_read_dict:   dictionary for additional motor read kwargs
        - ikwargs_dict:        dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}

        - mobj_dict:            dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:            same as mobj_dict but for instruments
        - vars:                 list of variable names to return and, if applicable, plot, must match something in measurd_header (see below). defaults to DEFAULT_VARS in configuration.py
        - metadata:             dictionary with any additional entries to add to metadata
        - xvar:                 variable to plot on x axis. defaults time
        - filename_head:        path to directory in which to save data. defaults ot current working directory of notebook
        - filename:             name of file
        - savefile:             if true, saves file
        - plot:                 if true, plots data in vars
        - close_devices:        if true, closes all devices in mobj_dict and iobj_dict
    '''

    # parse motor scanning information contained in sequence_list
    sequence_motors, mtargets, mkwargs_move, motors, tols = helper.capture_sequence_information(sequence_list)

    # initialize motors and instruments, and create kwarg dictionaries for each
    mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict = helper.initialize_instruments_and_motors(ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict)

    # create motor and instrument header
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
    motors_header_measured = [MOTOR_DICT[m]['name']+str(' measured') for m in ACTIVE_MOTORS] # for tracking measured motor positions
    header = ['Time (s)']+instruments_header+motors_header_measured

    # take default value for vars, and check that vars is valid
    if vars==[]:
        vars = DEFAULT_VARS
    if not set(vars).issubset(set(header)):
        raise ValueError(f'vars {vars} not in measurd values {header}')

    # setup file for writing - adds metadata at top and writes the data header
    if savefile:
        # setup metadata
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        metamotors_dict = dict([(MOTOR_DICT[m]['name'], measured_motors_dict[m]) for m in list(measured_motors_dict.keys())])
        metadata = {**metadata, **metamotors_dict}
        # setup and write metadta + header to file
        fname = helper.setup_filename(filename, filename_head, 'motor_sequence_'+'_'.join([f'{m}' for m in sequence_motors]))
        helper.write_file_header(fname, header, metadata)

    # setup plots for displaying
    if plot:
        xlabel = 'Time (s)'
        if xvar!='time':
            xlabel = MOTOR_DICT[xvar]['name']
        fig, axes, plot_handles_dict, xrange, vdata1d_dict = plotters.setup_1d_plots_append(vars, xlabel)

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    # loop
    t = 0
    tic = time.perf_counter()
    for ii, m in enumerate(sequence_motors):
        mtarget = mtargets[ii]
        move_dict = {m:(mtarget, False)}
        mkwargs_move_dict = {m:{**mkwargs_move[ii], 'check_stability':False}}
        #print(move_dict)
        #print(mobj_dict)
        #print(mkwargs_move_dict)
        #print(mkwargs_read_dict)
        mobj_dict = helper.move_motors(move_dict, mobj_dict, mkwargs_move_dict, mkwargs_read_dict, print_flag=False)
        tol = tols[ii]
        measured_motors_dict, mobj_dict = helper.read_motors([m], mobj_dict, mkwargs_read_dict)
        mcurrent = measured_motors_dict[m]
        #time.sleep(time_constant)
        while np.abs(mcurrent - mtarget) > tol:

            # if run is False, break inner loop
            if not run.locked_read():
                break

            # get new time
            toc = time.perf_counter()
            t = toc - tic

            # acquire data and read actual motor positions - should test how long this takes for typical motor setup
            #t0 = time.time()
            instrument_data_dict, iobj_dict = helper.read_instruments(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
            measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
            #tf = time.time()
            #print(tf-t0)

            # update mcurrent from measured_motors
            mcurrent = measured_motors_dict[m]

            # combine dictionaries to keep all new data in a single dictionary
            newdata = {**instrument_data_dict, **measured_motors_dict}

            # convert acquired data from dictionaries to list for file writing
            instrument_data = [instrument_data_dict[h] for h in instruments_header]
            measured_motors_data = [measured_motors_dict[m] for m in ACTIVE_MOTORS]

            # append new data to file
            if savefile:
                helper.append_data_to_file(fname, list([t])+instrument_data+measured_motors_data)

            # update plots
            if plot==True:
                newx = t
                if xvar!='time':
                    newx = newdata[xvar]
                xrange, vdata1d_dict = plotters.update_1d_plots_append(fig, axes, vars, plot_handles_dict, xrange, vdata1d_dict, newx, newdata)

        # if run is False, break outer loop
        if not run.locked_read():
            break

    # wait for thread to finish
    run.locked_update(False)
    user_input_thread.join()

    # close instruments and motors
    if close_devices:
        helper.close_instruments(iobj_dict)
        helper.close_motors(mobj_dict)
    else:
        return iobj_dict, mobj_dict

#####################
### Lab Utilities ###
#####################

def monitor_motors(motors, mkwargs_read_dict={}, filename_head=None, filename=None, savefile=False, plot=True):
    '''
    Utility to measure the value of a motor as function of time, with a quick way of exiting.

    args:
        - motors:                list of motors to display. single motor string will be turned into a list.
        - mkwargs_read_dict:     dictionary of motor read kwargs
        - filename_head:
        - filename:
        - savefile:
        - plot
    '''

    if type(motors)is str:
        motors = [motors]

    # setup file for writing
    if savefile:
        append=False
        if filename==None:
            append=True
        filename_head, filename = helper.generate_filename(filename_head, filename, 'monitor_motors')
        if append==True:
                for m in motors:
                    filename = filename+f'_{m}'
        fname = helper.get_unique_filename(filename_head, filename)
        header = ['Time (s)']+[MOTOR_DICT[m]['name'] for m in motors]
        helper.write_file_header(fname, header)

    # initialize motors
    mobj_dict = helper.initialize_motors(motors, {})
    mkwargs_read_dict = helper.construct_mkwargs_dict(mkwargs_read_dict)

    # setup plots for displaying
    if plot:
        xlabel = 'Time (s)'
        fig, axes, plot_handles_dict, xrange, vdata1d_dict = plotters.setup_1d_plots_append(motors, xlabel, figsize=(4,3))

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    t0 = time.time()
    while run.locked_read():

        measured_motors_dict, mobj_dict = helper.read_motors(motors, mobj_dict, mkwargs_read_dict)
        t = time.time() - t0

        measured_positions_data = list(measured_motors_dict.values())

        # add to file
        if savefile:
            helper.append_data_to_file(fname, [t]+measured_positions_data)

        if plot==True:
            xrange, vdata1d_dict = plotters.update_1d_plots_append(fig, axes, motors, plot_handles_dict, xrange, vdata1d_dict, t, measured_motors_dict)
        time.sleep(0.1)

    run.locked_update(False)
    user_input_thread.join()

    # close motors
    helper.close_motors(mobj_dict)

def list_motors():
    '''
    try to initialize motors in MOTOR_DICT and return those that do so properly. prints motors that fail to initialize.
    '''

    m_list = []
    for m in MOTOR_DICT.keys():
        try:
            mobj = MOTOR_DICT[m]['init']()
            m_list.append(m)
            MOTOR_DICT[m]['close'](mobj)
        except:
            print(f'motor {m} not able to initialize')

    return m_list

##############################
###                        ###
###   Optics Lab Methods   ###
###                        ###
##############################

#################################
### Rotations and Corotations ###
#################################

def rotate_scan(start_angle, end_angle, step_size, axis_index=1, mkwargs_read_dict={}, ikwargs_dict={}, mobj_dict={}, iobj_dict={}, vars=[], metadata={}, filename_head=None, filename=None, savefile=True, plot=True, close_devices=True):
    '''
    rotation scan for axes 1, 2 or 3

    args:
        - start_angle
        - end_angle
        - step_size
        - axis_index            axis to rotate. 1, 2, 3
        - mkwargs_read_dict:    dictionary of mkwargs for motor read functions
        - ikwargs_dict:            dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}}

        - mobj_dict:            dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:            same as mobj_dict but for instruments
        - vars:             list of variable names to return and, if applicable, plot, must match something in measurd_header (see below). defaults to DEFAULT_VARS in configuration.py
        - metadata:         dictionary with any additional entries to add to metadata
        - filename_head:        path to directory in which to save data. defaults ot current working directory of notebook
        - filename:             name of file
        - savefile:             if true, saves file
        - plot:                 if true, plots data in vars
        - close_devices:        if true, closes all devices in mobj_dict and iobj_dict
    '''

    # initialize motors and instruments, and create kwarg dictionaries for each
    mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict = helper.initialize_instruments_and_motors(ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict)

    # create motor and instrument header
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
    motors_header_measured = [MOTOR_DICT[m]['name']+str(' measured') for m in ACTIVE_MOTORS] # for tracking measured motor positions
    header = [f'Angle {axis_index} (deg)']+instruments_header+motors_header_measured

    # take default value for vars, and check that vars is valid
    if vars==[]:
        vars = DEFAULT_VARS
    if not set(vars).issubset(set(header)):
        raise ValueError(f'vars {vars} not in measurd values {header}')

    # setup file for writing - adds metadata at top and writes the data header
    if savefile:
        # setup metadata
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        metamotors_dict = dict([(MOTOR_DICT[m]['name'], measured_motors_dict[m]) for m in list(measured_motors_dict.keys())])
        metadata = {**metadata, **metamotors_dict}
        # setup and write metadta + header to file
        fname = helper.setup_filename(filename, filename_head, 'rotate_scan')
        helper.write_file_header(fname, header, metadata)

    # setup plots for displaying
    if plot:
        xlabel = f'Angle {axis_index} (deg)'
        fig, axes, plot_handles_dict, xrange, vdata1d_dict = plotters.setup_1d_plots_append(vars, xlabel)

    # select axis to rotate and setup both that axis and conjugate axis
    if axis_index not in [1,2,3]:
        raise ValueError('Invalid axis_index, please select either 1, 2, or 3')
    axis_name = f'axis_{axis_index}'
    move_axis = MOTOR_DICT[axis_name]['move']
    move_back = MOTOR_DICT[axis_name]['move_back']
    read_axis = MOTOR_DICT[axis_name]['read']
    axis_obj = mobj_dict[axis_name]

    # convert input to angle lists
    angles = helper.get_motor_range(start_angle, end_angle, step_size)

    # move axis to start
    axis_obj = move_axis(start_angle-move_back, axis=axis_obj)
    axis_obj = move_axis(start_angle, axis=axis_obj)
    mobj_dict[axis_name] = axis_obj

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    # scan
    for ii, angle in enumerate(angles):

        # if run is False, break loop
        if not run.locked_read():
            break

        # move
        axis_obj = move_axis(angle, axis=axis_obj)
        mobj_dict[axis_name] = axis_obj

        # acquire data and read actual motor positions - should test how long this takes for typical motor setup
        #t0 = time.time()
        instrument_data_dict, iobj_dict = helper.read_instruments(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        #tf = time.time()
        #print(tf-t0)

        # combine dictionaries to keep all new data in a single dictionary
        newdata = {**instrument_data_dict, **measured_motors_dict}

        # convert acquired data from dictionaries to list for file writing
        instrument_data = [instrument_data_dict[h] for h in instruments_header]
        measured_motors_data = [measured_motors_dict[m] for m in ACTIVE_MOTORS]

        # append new data to file
        if savefile:
            helper.append_data_to_file(fname, [angle]+instrument_data+measured_motors_data)

        # update plots
        if plot==True:
            xrange, vdata1d_dict = plotters.update_1d_plots_append(fig, axes, vars, plot_handles_dict, xrange, vdata1d_dict, angle, newdata)

    axis_obj = move_axis(angle, axis=axis_obj)
    mobj_dict[axis_name] = axis_obj

    # wait for thread to finish
    run.locked_update(False)
    user_input_thread.join()

    # close instruments and motors
    if close_devices:
        helper.close_instruments(iobj_dict)
        helper.close_motors(mobj_dict)
    else:
        return iobj_dict, mobj_dict

def corotate_scan(start_angle, end_angle, step_size, angle_offset, rate_axis_2=1, mkwargs_read_dict={}, ikwargs_dict={}, mobj_dict={}, iobj_dict={}, vars=[], metadata={}, filename_head=None, filename=None, savefile=True, plot=True, close_devices=True):
    '''
    Corotation scan moving axes 1 and 2, typically representing wave plates.

    axis_2_angle = rate_axis_2*axis_1_angle + angle_offset

    args:
        - start_angle
        - end_angle
        - step_size
        - angle_offset      offset angle between
        - rate_axis_2       rate at which axis_2 moves relative to axis_1. can be negative.
        - mkwargs_read_dict:    dictionary of mkwargs for motor read functions
        - ikwargs_dict:        dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}}

        - mobj_dict:        dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:        same as mobj_dict but for instruments
        - vars:             list of variable names to return and, if applicable, plot, must match something in measurd_header (see below). defaults to DEFAULT_VARS in configuration.py
        - metadata:         dictionary with any additional entries to add to metadata
        - filename_head:    path to directory in which to save data. defaults ot current working directory of notebook
        - filename:         name of file
        - savefile:         if true, saves file
        - plot:             if true, plots data in vars
        - close_devices:    if true, closes all devices in mobj_dict and iobj_dict
    '''

    # initialize motors and instruments, and create kwarg dictionaries for each
    mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict = helper.initialize_instruments_and_motors(ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict)

    # create motor and instrument header
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
    motors_header_measured = [MOTOR_DICT[m]['name']+str(' measured') for m in ACTIVE_MOTORS] # for tracking measured motor positions
    header = [f'Angle 1 (deg)', 'Angle 2 (deg)']+instruments_header+motors_header_measured

    # take default value for vars, and check that vars is valid
    if vars==[]:
        vars = DEFAULT_VARS
    if not set(vars).issubset(set(header)):
        raise ValueError(f'vars {vars} not in measurd values {header}')

    # setup file for writing - adds metadata at top and writes the data header
    if savefile:
        # setup metadata
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        metamotors_dict = dict([(MOTOR_DICT[m]['name'], measured_motors_dict[m]) for m in list(measured_motors_dict.keys())])
        metadata = {**metadata, **metamotors_dict}
        # setup and write metadta + header to file
        fname = helper.setup_filename(filename, filename_head, 'corotate_scan')
        helper.write_file_header(fname, header, metadata)


    # setup plots for displaying
    if plot:
        xlabel = f'Angle 1 (deg)'
        fig, axes, plot_handles_dict, xrange, vdata1d_dict = plotters.setup_1d_plots_append(vars, xlabel)

    # select axis to rotate and setup both that axis and conjugate axis
    axis_1 = mobj_dict['axis_1']
    axis_2 = mobj_dict['axis_2']
    move_axis_1 = MOTOR_DICT['axis_1']['move']
    move_axis_2 = MOTOR_DICT['axis_2']['move']
    read_axis_1 = MOTOR_DICT['axis_1']['read']
    read_axis_2 = MOTOR_DICT['axis_2']['read']
    move_back_1 = MOTOR_DICT['axis_1']['move_back']
    move_back_2 = MOTOR_DICT['axis_2']['move_back']

    # convert input to angle lists
    angles_1 = helper.get_motor_range(start_angle, end_angle, step_size)
    angles_2 = rate_axis_2*angles_1 + angle_offset

    # move axis to start
    axis_1, axis_2 = newport.corotate_axes(1, 2, start_angle-move_back_1, start_angle+angle_offset-move_back_2, axis_1=axis_1, axis_2=axis_2)
    axis_1, axis_2 = newport.corotate_axes(1, 2, start_angle, start_angle+angle_offset, axis_1=axis_1, axis_2=axis_2)
    mobj_dict['axis_1'] = axis_1
    mobj_dict['axis_2'] = axis_2

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    # scan
    for ii, angle in enumerate(angles_1):

        # if run is False, break loop
        if not run.locked_read():
            break

        # move angles
        angle_1 = angle
        angle_2 = angles_2[ii]
        axis_1, axis_2 = newport.corotate_axes(1, 2, angle, angle+angle_offset, axis_1=axis_1, axis_2=axis_2)
        mobj_dict['axis_1'] = axis_1
        mobj_dict['axis_2'] = axis_2

        # acquire data and read actual motor positions - should test how long this takes for typical motor setup
        #t0 = time.time()
        instrument_data_dict, iobj_dict = helper.read_instruments(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
        measured_motors_dict, mobj_dict = helper.read_motors(ACTIVE_MOTORS, mobj_dict, mkwargs_read_dict)
        #tf = time.time()
        #print(tf-t0)

        # combine dictionaries to keep all new data in a single dictionary
        newdata = {**instrument_data_dict, **measured_motors_dict}

        # convert acquired data from dictionaries to list for file writing
        instrument_data = [instrument_data_dict[h] for h in instruments_header]
        measured_motors_data = [measured_motors_dict[m] for m in ACTIVE_MOTORS]

        # append new data to file
        if savefile:
            helper.append_data_to_file(fname, [angle_1, angle_2]+instrument_data+measured_motors_data)

        # update plots
        if plot==True:
            xrange, vdata1d_dict = plotters.update_1d_plots_append(fig, axes, vars, plot_handles_dict, xrange, vdata1d_dict, angle, newdata)

    # move motor back to original positions
    axis_1, axis_2 = newport.corotate_axes(1, 2, start_angle-move_back_1, start_angle+angle_offset-move_back_2, axis_1=axis_1, axis_2=axis_2)
    mobj_dict['axis_1'] = axis_1
    mobj_dict['axis_2'] = axis_2
    axis_1, axis_2 = newport.corotate_axes(1, 2, start_angle, start_angle+angle_offset, axis_1=axis_1, axis_2=axis_2)
    mobj_dict['axis_1'] = axis_1
    mobj_dict['axis_2'] = axis_2

    # wait for thread to finish
    run.locked_update(False)
    user_input_thread.join()

    # close instruments and motors
    if close_devices:
        helper.close_instruments(iobj_dict)
        helper.close_motors(mobj_dict)
    else:
        return iobj_dict, mobj_dict

def rotate_map(map_dict, start_angle, end_angle, step_size, mkwargs_read_dict={}, ikwargs_dict={}, mobj_dict={}, iobj_dict={}, vars=[], metadata={}, filename_head=None, filename=None, savefile=True, print_flag=False, plot=False, close_devices=True):
    '''
    method to record single axis rotation measurements on instruments as a function of motors specified by dictionary map_dict. meant to acquire data one point at a time in a multidimensional motor space.

    rotation scan for axes 1, 2 or 3

    args:
        - map_dict:         dictionary of key:value pairs where key is name of a motor in MOTOR_DICT and value is a tuple (start, stop step, mkwargs_move) or (positions, mkwargs_move), and where mkwargs_move is a dictionary where keys are kwarg names for the motor move function, and value the corresponding value to set during move calls. motors are scanned from left to right.

            ex: {'temp':(100,200,10,{'wait_time':30})} or {'temp':(np.arange(100,200,10),{'wait_time':30})}

        - start_angle:
        - end_angle:
        - step_size:
        - mkwargs_read_dict:    dictionary of mkwargs for motor read functions. if left empty, defaults to no kwargs
        - ikwargs_dict:        dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}}

        - mobj_dict:        dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:        same as mobj_dict but for instruments
        - vars:             list of variable names to return and, if applicable, plot, must match something in measurd_header (see below). defaults to DEFAULT_VARS in configuration.py
        - metadata:         dictionary with any additional entries to add to metadata
        - filename_head:    path to directory in which to save data. defaults ot current working directory of notebook
        - filename:         name of file
        - savefile:         if true, saves file
        - print_flag:       if true, print motor positions during scan
        - plot:             if true, plots data in vars
        - close_devices:    if true, closes all devices in mobj_dict and iobj_dict
    '''

    # parse motor scanning information contained in map_dict
    motors, mranges, mkwargs_move_dict = helper.capture_motor_information(map_dict)

    # initialize motors and instruments, and create kwarg dictionaries for each
    mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict = helper.initialize_instruments_and_motors(ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict)
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)

    # take default value for vars, and check that vars is valid
    if vars==[]:
        vars = DEFAULT_VARS
    measured_vars = [MOTOR_DICT[m]['name'] for m in ACTIVE_MOTORS]+instruments_header
    if not set(vars).issubset(set(measured_vars)):
        raise ValueError(f'vars {vars} not in measurd values {measured_vars}')

    # setup filenames
    filename = helper.setup_filename(filename, filename_head, 'rotate_map')

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    # generate motor scan positions recursively.
    # positions: list of positions where each element contains positions of each motor for nth step in scan
    positions = helper.gen_positions_recurse(mranges, len(mranges)-1)
    num_pos = len(positions)

    # move motors to start
    move_dict = {}
    for jj, m in enumerate(motors):
        move_dict[m] = (positions[0][jj], True)
    mobj_dict = helper.move_motors(move_dict, mobj_dict, mkwargs_move_dict, mkwargs_read_dict, print_flag=print_flag)

    # loop over positions, only moving a motor if its target position has changed.
    start_pos = positions[0]
    current_pos = start_pos
    for ii in tqdm(range(num_pos)):
        pos = positions[ii]

        # if run is False, break loop
        if not run.locked_read():
            break

        # move motors for which position has changed, moving first to "go back" position if moving back to start of initial motor position
        move_dict = {}
        for jj, m in enumerate(motors):
            mtarget = pos[jj]
            if current_pos[jj]!=mtarget:
                move_back_flag=False
                if start_pos[jj]==mtarget:
                    move_back_flag=True
                move_dict[m] = (mtarget, move_back_flag)
        mobj_dict = helper.move_motors(move_dict, mobj_dict, mkwargs_move_dict, mkwargs_read_dict, print_flag=print_flag)
        current_pos = pos

        # setup each filename
        expanded_filename = filename
        for jj, m in enumerate(motors):
            p = pos[jj]
            expanded_filename = expanded_filename+f'_{m}{p}'

        # scan
        iobj_dict, mobj_dict = rotate_scan(start_angle, end_angle, step_size, axis_index=1, mkwargs_read_dict=mkwargs_read_dict, ikwargs_dict=ikwargs_dict, mobj_dict=mobj_dict, iobj_dict=iobj_dict, vars=vars, metadata=metadata, filename_head=filename_head, filename=filename, savefile=True, plot=plot, close_devices=False)

    # wait for thread to finish
    run.locked_update(False)
    user_input_thread.join()

    # close instruments and motors
    if close_devices:
        helper.close_instruments(iobj_dict)
        helper.close_motors(mobj_dict)
    else:
        return iobj_dict, mobj_dict

def corotate_map(map_dict, start_angle, end_angle, step_size, angle_offset, rate_axis_2=1, mkwargs_read_dict={}, ikwargs_dict={}, mobj_dict={}, iobj_dict={}, vars=[], metadata={}, filename_head=None, filename=None, savefile=True, print_flag=False, plot=False, close_devices=True):
    '''
    method to record corotation measurements on instruments as a function of motors specified by dictionary map_dict. meant to acquire data one point at a time in a multidimensional motor space.

    Corotation scan moving axes 1 and 2, typically representing wave plates.

    axis_2_angle = rate_axis_2*axis_1_angle + angle_offset

    args:
        - map_dict:         dictionary of key:value pairs where key is name of a motor in MOTOR_DICT and value is a tuple (start, stop step, mkwargs_move) or (positions, mkwargs_move), and where mkwargs_move is a dictionary where keys are kwarg names for the motor move function, and value the corresponding value to set during move calls. motors are scanned from left to right.

            ex: {'temp':(100,200,10,{'wait_time':30})} or {'temp':(np.arange(100,200,10),{'wait_time':30})}

        - start_angle:
        - end_angle:
        - step_size:
        - angle_offset:
        - mkwargs_read_dict:    dictionary of mkwargs for motor read functions. if left empty, defaults to no kwargs
        - ikwargs_dict:        dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}}

        - mobj_dict:        dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:        same as mobj_dict but for instruments
        - vars:             list of variable names to return and, if applicable, plot, must match something in measurd_header (see below). defaults to DEFAULT_VARS in configuration.py
        - metadata:         dictionary with any additional entries to add to metadata
        - filename_head:    path to directory in which to save data. defaults ot current working directory of notebook
        - filename:         name of file
        - savefile:         if true, saves file
        - print_flag:       if true, print motor positions during scan
        - plot:             if true, plots data in vars
        - close_devices:    if true, closes all devices in mobj_dict and iobj_dict
    '''

    # parse motor scanning information contained in map_dict
    motors, mranges, mkwargs_move_dict = helper.capture_motor_information(map_dict)

    # initialize motors and instruments, and create kwarg dictionaries for each
    mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict = helper.initialize_instruments_and_motors(ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, mkwargs_read_dict, ikwargs_dict, mobj_dict, iobj_dict)
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)

    # take default value for vars, and check that vars is valid
    if vars==[]:
        vars = DEFAULT_VARS
    measured_vars = [MOTOR_DICT[m]['name'] for m in ACTIVE_MOTORS]+instruments_header
    if not set(vars).issubset(set(measured_vars)):
        raise ValueError(f'vars {vars} not in measurd values {measured_vars}')

    # setup filenames
    filename = helper.setup_filename(filename, filename_head, 'corotate_map')

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    # generate motor scan positions recursively.
    # positions: list of positions where each element contains positions of each motor for nth step in scan
    positions = helper.gen_positions_recurse(mranges, len(mranges)-1)
    num_pos = len(positions)

    # move motors to start
    move_dict = {}
    for jj, m in enumerate(motors):
        move_dict[m] = (positions[0][jj], True)
    mobj_dict = helper.move_motors(move_dict, mobj_dict, mkwargs_move_dict, mkwargs_read_dict, print_flag=print_flag)

    # loop over positions, only moving a motor if its target position has changed.
    start_pos = positions[0]
    current_pos = start_pos
    for ii in tqdm(range(num_pos)):
        pos = positions[ii]

        # if run is False, break loop
        if not run.locked_read():
            break

        # move motors for which position has changed, moving first to "go back" position if moving back to start of initial motor position
        move_dict = {}
        for jj, m in enumerate(motors):
            mtarget = pos[jj]
            if current_pos[jj]!=mtarget:
                move_back_flag=False
                if start_pos[jj]==mtarget:
                    move_back_flag=True
                move_dict[m] = (mtarget, move_back_flag)
        mobj_dict = helper.move_motors(move_dict, mobj_dict, mkwargs_move_dict, mkwargs_read_dict, print_flag=print_flag)
        current_pos = pos

        # setup each filename
        expanded_filename = filename
        for jj, m in enumerate(motors):
            p = pos[jj]
            expanded_filename = expanded_filename+f'_{m}{p}'

        # scan
        iobj_dict, mobj_dict = corotate_scan(start_angle, end_angle, step_size, angle_offset, rate_axis_2=rate_axis_2, mkwargs_read_dict=mkwargs_read_dict, ikwargs_dict=ikwargs_dict, mobj_dict=mobj_dict, iobj_dict=iobj_dict, vars=vars, metadata=metadata, filename_head=filename_head, filename=filename, savefile=True, plot=plot, close_devices=False)

    # wait for thread to finish
    run.locked_update(False)
    user_input_thread.join()

    # close instruments and motors
    if close_devices:
        helper.close_instruments(iobj_dict)
        helper.close_motors(mobj_dict)
    else:
        return iobj_dict, mobj_dict

#########################
### Balancing Methods ###
#########################

def find_balance_angle(start_angle, end_angle, step_size, balance_at=0, offset=0, bal_axis=2, go_to_balance_angle=True, var='Demod 1 x', ikwargs_dict={}, mobj_dict={}, iobj_dict={}):
    '''
    Assuming we are measuring in DC mode above a transition or on GaAs, carries out a rotate_scan. Find angle by carrying out a linear fit, such that the angle range should be taken to be very small.

    By default, moves axis 2 and keep axis 1 static to find balance angle.

    If go_to_balance_angle is set to true, moves axis to balance angle after scan

    args:
        - start_angle
        - end_angle
        - step_size
        - balance_at:           angle of static axis at which to balance
        - offset:               target balanced signal
        - bal_axis:             axis with which to balance
        - go_to_balance_angle:  if true, set axis to balance angle
        - var:                  variable signal to monitor
        - ikwargs_dict:            dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}}

        - mobj_dict:            dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:            same as mobj_dict but for instruments

    return:
        - balance_angle:    approximate balance angle from fit
        - slope:            slope of signal vs angle, to be used in autobalancing method

    '''

    if bal_axis==2:
        bal_axis='axis_2'
        static_axis='axis_1'
    else:
        bal_axis='axis_1'
        static_axis='axis_2'

    # capture instrument information and initialize
    # ikwargs_dict: dictionary of kwargs for each instrument in ACTIVE_INSTRUMENTS, including kwargs specified in ikwargs_dict. if ikwargs_dict is empty, defaults to {'inst1':{},....,'instN':{}} for N instruments
    # iobjs_dict: dictionary of instrument handle objects for each instrument in instruments
    # instruments_header: list of strings with instrument measurement names for each instrument.
    ikwargs_dict = helper.construct_ikwargs_dict(ikwargs_dict)
    iobj_dict = helper.initialize_instruments(ACTIVE_INSTRUMENTS, iobj_dict)
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)

    # initialize ACTIVE_MOTORS to read during measurement
    # mobj_dict: dictionary of motor handle objects for each motor in motors
    mobj_dict = helper.initialize_motors(['axis_1', 'axis_2'], mobj_dict)

    # define axes objects and funcs
    static_axis_obj = mobj_dict[static_axis]
    bal_axis_obj = mobj_dict[bal_axis]
    move_static_axis = MOTOR_DICT[static_axis]['move']
    move_bal_axis = MOTOR_DICT[bal_axis]['move']
    move_back_static_axis = MOTOR_DICT[static_axis]['move_back']
    move_back_bal_axis = MOTOR_DICT[bal_axis]['move_back']

    # move both motors to balance_at
    static_axis_obj = move_static_axis(balance_at-move_back_static_axis, static_axis_obj)
    bal_axis_obj = move_bal_axis(balance_at-move_back_bal_axis, bal_axis_obj)
    static_axis_obj = move_static_axis(balance_at, static_axis_obj)
    bal_axis_obj = move_bal_axis(balance_at, bal_axis_obj)

    # convert input to angle lists
    angles = helper.get_motor_range(balance_at+start_angle, balance_at+end_angle, step_size)

    # rotate axis and collect data
    signal = np.zeros(len(angles))

    # setup figure for plotting
    xlabel = f'Angle (deg)'
    fig, axes, plot_handles_dict, xrange, vdata1d_dict = plotters.setup_1d_plots_append([var], xlabel)

    for ii, angle in enumerate(angles):

        # move
        bal_axis_obj = move_bal_axis(angle, bal_axis_obj)

        instrument_data_dict, iobj_dict = helper.read_instruments(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)

        signal[ii] = instrument_data_dict[var]

        # update plots
        xrange, vdata1d_dict = plotters.update_1d_plots_append(fig, axes, [var], plot_handles_dict, xrange, vdata1d_dict, angle, instrument_data_dict)

    # linear fit
    fitf = lambda x, m, b: m*x+b
    popt, pcov = opt.curve_fit(fitf, angles, signal-offset)
    balance_angle = -popt[1]/popt[0]
    slope = popt[0]
    angles_vect = np.linspace(start_angle+balance_at, balance_at+end_angle, 1000)
    fit = popt[0]*angles_vect + popt[1]

    # display fit result
    axes[0].plot(angles_vect, fit, '-', color='black')

    # set balance angle relative to 0
    balance_angle = balance_angle - balance_at

    if go_to_balance_angle:
        bal_axis_obj = move_bal_axis(balance_angle+balance_at, bal_axis_obj)

    mobj_dict[static_axis] = static_axis_obj
    mobj_dict[bal_axis] = bal_axis_obj
    helper.close_motors(mobj_dict)

    print(f'Balance angle: {balance_angle}')
    return balance_angle, slope

def autobalance(slope, tolerance, offset=0, bal_axis=2, var='Demod 1 x', ikwargs_dict={}, mobj_dict={}, iobj_dict={}, print_flag=True, close_devices=True):

    '''
    balance by minimizing signal (var-offset) below tolerance tol. By default, moves axis 2 and keep axis 1 static to find balance angle.

    args:
        - slope:                slope of var signal near balance angle, extracted from find_balance_angle for example
        - tol erance            tolerance for minizing signal. when signal below tol threshold, functions return.
        - offset:               target balanced signal
        - bal_axis:             axis with which to balance
        - var:                  variable signal to monitor
        - ikwargs_dict:            dictionary of key value pairs where keys are name of instruments in INSTRUMENT_DICT and values are dictionary of kwargs for the instrument read function. Defaults to instruments in ACTIVE_INSTRUMENTS with no kwargs.

            ex: {'zurich_lockin':{'time_constant':0.3, 'channel_index':3}}

        - mobj_dict:            dictionary containing key:value pairs where keys are motor name and value are motor handle objects
        - iobj_dict:            same as mobj_dict but for instruments
        - print_flag:           if true, print out balance angle

    return:
        - balance_angle:        angle which minizes

    '''

    if bal_axis==2:
        bal_axis='axis_2'
        static_axis='axis_1'
    else:
        bal_axis='axis_1'
        static_axis='axis_2'

    # capture instrument information and initialize
    # ikwargs_dict: dictionary of kwargs for each instrument in ACTIVE_INSTRUMENTS, including kwargs specified in ikwargs_dict. if ikwargs_dict is empty, defaults to {'inst1':{},....,'instN':{}} for N instruments
    # iobjs_dict: dictionary of instrument handle objects for each instrument in instruments
    # instruments_header: list of strings with instrument measurement names for each instrument.
    ikwargs_dict = helper.construct_ikwargs_dict(ikwargs_dict)
    iobj_dict = helper.initialize_instruments(ACTIVE_INSTRUMENTS, iobj_dict)
    instruments_header, iobj_dict = helper.get_instruments_header(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)

    # initialize ACTIVE_MOTORS to read during measurement
    # mobj_dict: dictionary of motor handle objects for each motor in motors
    mobj_dict = helper.initialize_motors(['axis_1', 'axis_2'], mobj_dict)

    # define axes objects and funcs
    static_axis_obj = mobj_dict[static_axis]
    bal_axis_obj = mobj_dict[bal_axis]
    read_static_axis = MOTOR_DICT[static_axis]['read']
    read_bal_axis = MOTOR_DICT[bal_axis]['read']
    move_static_axis = MOTOR_DICT[static_axis]['move']
    move_bal_axis = MOTOR_DICT[bal_axis]['move']
    move_back_static_axis = MOTOR_DICT[static_axis]['move_back']
    move_back_bal_axis = MOTOR_DICT[bal_axis]['move_back']

    # obtain the angle of static axis
    balance_at, bal_static_obj = read_static_axis(static_axis_obj)

    # minimize var-offset
    signal = 10000
    curr_pos, bal_axis_obj = read_bal_axis(bal_axis_obj)
    while (np.abs(signal-offset)>tolerance):


        data_dict, iobj_dict = helper.read_instruments(ACTIVE_INSTRUMENTS, iobj_dict, ikwargs_dict)
        signal = data_dict[var]
        new_pos = (curr_pos-(1/slope)*(signal-offset))%360
        bal_axis_obj = move_bal_axis(new_pos, bal_axis_obj)
        curr_pos = new_pos

    if print_flag:
        print(f'Balanced PID at {balance_at}. Balance angle: {curr_pos-balance_at}.')

    mobj_dict[static_axis] = static_axis_obj
    mobj_dict[bal_axis] = bal_axis_obj
    mobj_dict[static_axis] = static_axis_obj

    if close_devices:
        helper.close_motors(mobj_dict)

    return curr_pos-balance_at

#######################
### Other Utilities ###
#######################

def align_delay_stage(wait_time=5, range=(-125,125)):
    '''
    helpful tool for aligning 'delay stage' motor by continuously moving it from one point of the range to the other.

    args:
        - wait_time:        time to wait at each end point
        - range:            positions of delay stage to move between.
    '''

    # initialize delay stage:
    delay_stage = MOTOR_DICT['delay_stage']['init']()
    move = MOTOR_DICT['delay_stage']['move']
    close = MOTOR_DICT['delay_stage']['close']

    # start user input thread
    run = LockedVar(True)
    user_input_thread = threading.Thread(target=concurrency.press_esc_to_stop, args=(run,))
    user_input_thread.start()

    while run.locked_read():
        delay_stage = move(range[0], delay_stage)
        time.sleep(wait)
        delay_stage = move(range[1], delay_stage)
        time.sleep(wait)

    # wait for thread to finish
    run.locked_update(False)
    user_input_thread.join()

    # close delay stage
    close(delay_stage)
