'''
Main file for controlling lab equipment and orchestrating measurements, with a specific eye to procedures
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
import OrensteinLab_git.Measurement.Alex.control.newport as newport
import OrensteinLab_git.Measurement.Alex.control.lakeshore as ls
import OrensteinLab_git.Measurement.Alex.control.zurich as zurich
import OrensteinLab_git.Measurement.Alex.control.opticool as oc
import OrensteinLab_git.Measurement.Alex.control.razorbill as razorbill
import OrensteinLab_git.Measurement.Alex.control.attocube as atto
from OrensteinLab_git.Measurement.Alex.control.concurrency_classes import LockedVar, StoppableThread
from OrensteinLab_git.Measurement.Alex.motors import motor_dict, instrument_dict, meta_motors
from orenstein_analysis.measurement import loader, process

'''
Features to add:
    - overhaul how I handle lockin data acquisition and possible link this to metadata
'''

lockin_header = ['Demod x', 'Demod y', 'r', 'Demod x_R', 'Demod y_R', 'Demod r_R'] # this should be moved to zurich file somehow, ie, the output of reading zurich itself should tell us what the header are. Also I should make it so that all channels are saved and channel index just determines what gets plotted.

################
### Medthods ###
################

def lockin_time_series(recording_time, filename_head=None, filename=None, mobj_measure_dict={}, metadata={}, time_constant=0.3, channel_index=1, savefile=True, daq_objs=None, plot_flag=True, override_metadata=False):
    '''
    aquires data on the lockin over a specified length of time.
    '''

    # initialize zurich lockin and setup read function
    if daq_objs==None:
        daq_objs  = instrument_dict['zurich_lockin']['init']()
    read_lockin = instrument_dict['zurich_lockin']['read']
    lockin_header = list(read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index).keys())

    # initialize additional motors to measure during scan
    passed_measure_motors = True
    if mobj_measure_dict=={}:
        passed_measure_motors = False
        mobj_measure_dict = initialize_motors(meta_motors)
    measure_motors = list(mobj_measure_dict.keys())

    # setup metadata - ie, for quick reference of starting state before measurement
    if override_metadata==True:
        pass
    else:
        metadata = {**metadata, **generate_metadata(mobj_measure_dict)}

    # initialize data bins
    time_record = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # setup plot
    if plot_flag==True:
        fig, axes = plt.subplots(3,1, figsize=(8,10))
        y_labels = ['Demod x', 'Demod y', 'R']
        for ii, ax in enumerate(axes):
            ax.set_xlabel('Time (s)')
            ax.set_ylabel(y_labels[ii])
            ax.grid(True)
        draw_x, = axes[0].plot([],'-o')
        draw_y, = axes[1].plot([],'-o')
        draw_r, = axes[2].plot([],'-o')
        fig.canvas.draw()
        fig.show()

    # setup file for writing
    if savefile:
        filename_head, filename = generate_filename(filename_head, filename, 'lockin_time_series')
        fname = get_unique_filename(filename_head, filename)
        header_motors = [motor_dict[m]['name'] for m in measure_motors]
        header = ['Time (s)']+lockin_header+header_motors
        write_file_header(fname, header, metadata)

    # loop
    t_delay = 0
    tic = time.perf_counter()
    while (t_delay<recording_time):
        time.sleep(time_constant)
        toc = time.perf_counter()
        t_delay = toc - tic

        lockin_data = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index)
        x = lockin_data['Demod x']
        y = lockin_data['Demod y']
        r = lockin_data['Demod r']

        # read additional motors
        measured_positions_dict = read_motors(mobj_measure_dict)
        measured_positions = [measured_positions_dict[m] for m in measure_motors]

        time_record = np.append(time_record, t_delay)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # update plot
        if plot_flag==True:
            draw_x.set_data(time_record-time_record[0],demod_x)
            draw_y.set_data(time_record-time_record[0],demod_y)
            draw_r.set_data(time_record-time_record[0],demod_r)
            for ax in axes:
                ax.relim()
                ax.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

        # update file
        if savefile:
            vars = [t_delay]+list(lockin_data.values())+measured_positions
            append_data_to_file(fname, vars)

    # close motors
    if not passed_measure_motors:
        close_motors(mobj_measure_dict)

    return time_record, demod_x, demod_y, demod_r

def rotate_scan(start_angle, end_angle, step_size, filename_head=None, filename=None, axis_index=1, mobj_measure_dict={}, metadata={}, showplot=True, time_constant=0.3, channel_index=1, R_channel_index=4, daq_objs=None, axis_1=None, axis_2=None, savefile=True, override_metadata=False):

    # initialize zurich lockin and setup read function
    if daq_objs==None:
        init_func = instrument_dict['zurich_lockin']['init']
        daq_objs = init_func()
    read_lockin = instrument_dict['zurich_lockin']['read']
    lockin_header = list(read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index).keys())

    # initialize axes and setup move functions
    if axis_1 == None:
        init_func = motor_dict['axis_1']['init']
        axis_1 = init_func()
    if axis_2 == None:
        init_func = motor_dict['axis_2']['init']
        axis_2 = init_func()
    move_axis_1 = motor_dict['axis_1']['move']
    move_axis_2 = motor_dict['axis_2']['move']
    read_axis_1 = motor_dict['axis_1']['read']
    read_axis_2 = motor_dict['axis_2']['read']
    move_back_1 = motor_dict['axis_1']['move_back']
    move_back_2 = motor_dict['axis_2']['move_back']
    if axis_index==1:
        axis = axis_1
        move_axis = motor_dict['axis_1']['move']
        read_axis = motor_dict['axis_1']['read']
        move_back = motor_dict['axis_1']['move_back']
        other_axis = axis_2
        move_other_axis = motor_dict['axis_2']['move']
        read_other_axis = motor_dict['axis_2']['read']
        move_other_back = motor_dict['axis_2']['move_back']
    elif axis_index==2:
        axis = axis_2
        move_axis = motor_dict['axis_2']['move']
        read_axis = motor_dict['axis_2']['read']
        move_back = motor_dict['axis_2']['move_back']
        other_axis = axis_1
        move_other_axis = motor_dict['axis_1']['move']
        read_other_axis = motor_dict['axis_1']['read']
        move_other_back = motor_dict['axis_1']['move_back']
    else:
        raise ValueError('Invalid axis_index, please select either 1 or 2.')

    # initialize additional motors to measure during scan
    passed_measure_motors = True
    if mobj_measure_dict=={}:
        passed_measure_motors = False
        mobj_measure_dict = initialize_motors(list(set(meta_motors) - set(['axis_1', 'axis_2'])))
    measure_motors = list(mobj_measure_dict.keys())

    # setup metadata - ie, for quick reference of starting state before measurement
    if override_metadata==True:
        pass
    else:
        metadata = {**metadata, **generate_metadata(mobj_measure_dict)}

    # setup measureables
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # convert input to angle lists
    angles = get_motor_range(start_angle, end_angle, step_size)

    # setup file for writing
    if savefile:
        filename_head, filename = generate_filename(filename_head, filename, 'rotate_scan')
        fname = get_unique_filename(filename_head, filename)
        header_motors = [motor_dict[m]['name']+str(' measured') for m in measure_motors]
        header = [motor_dict['axis_1']['name'], motor_dict['axis_2']['name']]+lockin_header+header_motors
        write_file_header(fname, header, metadata)

    # setup plot
    if showplot==True:
        fig, axes = plt.subplots(3, 1, figsize=(8,10))
        y_labels = ['Demod x', 'Demod y', 'R']
        for ii, ax in enumerate(axes):
            ax.set_xlabel('Angle (deg)')
            ax.set_ylabel(y_labels[ii])
            ax.grid(True)
        draw_x, = axes[0].plot([],'-o')
        draw_y, = axes[1].plot([],'-o')
        draw_r, = axes[2].plot([],'-o')
        fig.canvas.draw()
        fig.show()

    # scan
    for ii, angle in enumerate(angles):
        if (angle == start_angle):
            move_axis(angle-move_back, axis=axis)
            #move_other_axis(angle-move_other_back, axis=other_axis)
            #move_other_axis(angle, axis=other_axis)
            move_axis(angle, axis=axis)
            time.sleep(2)
        move_axis(angle, axis=axis)
        time.sleep(time_constant)

        # read lockin and rotators
        lockin_data = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index)
        angle_pos_1 = read_axis_1(axis=axis_1, print_flag=False)
        angle_pos_2 = read_axis_2(axis=axis_2, print_flag=False)
        if axis_index==1:
            angle_pos = angle_pos_1
        elif axis_index==2:
            angle_pos = angle_pos_2

        # extract lockin data
        x = lockin_data['Demod x']
        y = lockin_data['Demod y']
        r = lockin_data['Demod r']

        # read additional motors
        measured_positions_dict = read_motors(mobj_measure_dict)
        measured_positions = [measured_positions_dict[m] for m in measure_motors]

        position = np.append(position, angle_pos)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # update plot
        if showplot == True:
            draw_x.set_data(position, demod_x)
            draw_y.set_data(position, demod_y)
            draw_r.set_data(position, demod_r)
            for ax in axes:
                ax.relim()
                ax.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

        # write to file
        if savefile:
            vars = [angle_pos_1, angle_pos_2]+list(lockin_data.values())+measured_positions
            append_data_to_file(fname, vars)

    # move motors back to original positions
    move_axis(start_angle, axis=axis)
    #move_other_axis(start_angle, axis=other_axis)

    # close motors
    if not passed_measure_motors:
        close_motors(mobj_measure_dict)

    return position, demod_x, demod_y, demod_r

def corotate_scan(start_angle, end_angle, step_size, angle_offset, rate_axis_2=1, filename_head=None, filename=None, mobj_measure_dict={}, metadata={}, showplot=True, time_constant=0.3, channel_index=1, R_channel_index=4, daq_objs=None, axis_1=None, axis_2=None, savefile=True, override_metadata=False):
    '''
    Takes a corotation scan moving axes 1 and 2, typically representing half wave plates. axis 2 can be specificied to move at a rate greater than axis 1, such that axis_2_move = rate*axis_1_move

    To do:
        - build in ability to change scan direction

    '''

    # initialize zurich lockin and setup read function
    if daq_objs==None:
        init_func = instrument_dict['zurich_lockin']['init']
        daq_objs = init_func()
    read_lockin = instrument_dict['zurich_lockin']['read']
    lockin_header = list(read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index).keys())

    # initialize axes and setup move functions
    if axis_1 == None:
        init_func = motor_dict['axis_1']['init']
        axis_1 = init_func()
    if axis_2 == None:
        init_func = motor_dict['axis_2']['init']
        axis_2 = init_func()
    move_axis_1 = motor_dict['axis_1']['move']
    move_axis_2 = motor_dict['axis_2']['move']
    read_axis_1 = motor_dict['axis_1']['read']
    read_axis_2 = motor_dict['axis_2']['read']
    move_back_1 = motor_dict['axis_1']['move_back']
    move_back_2 = motor_dict['axis_2']['move_back']

    # initialize additional motors to measure during scan
    passed_measure_motors = True
    if mobj_measure_dict=={}:
        passed_measure_motors = False
        mobj_measure_dict = initialize_motors(list(set(meta_motors) - set(['axis_1', 'axis_2'])))
    measure_motors = list(mobj_measure_dict.keys())

    # setup metadata - ie, for quick reference of starting state before measurement
    if override_metadata==True:
        pass
    else:
        metadata = {**metadata, **generate_metadata(mobj_measure_dict)}

    # setup measureables
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # convert input to angle lists
    angles_1 = get_motor_range(start_angle, end_angle, step_size)
    angles_2 = rate_axis_2*angles_1 + angle_offset

    # setup file for writing
    if savefile:
        filename_head, filename = generate_filename(filename_head, filename, 'corotate_scan')
        fname = get_unique_filename(filename_head, filename)
        header_motors = [motor_dict[m]['name']+str(' measured') for m in measure_motors]
        header = header = [motor_dict['axis_1']['name'], motor_dict['axis_2']['name']]+lockin_header+header_motors
        write_file_header(fname, header, metadata)

    # setup plot
    if showplot==True:
        fig, axes = plt.subplots(3, 1, figsize=(8,10))
        y_labels = ['Demod x', 'Demod y', 'R']
        for ii, ax in enumerate(axes):
            ax.set_xlabel('Angle_1 (deg)')
            ax.set_ylabel(y_labels[ii])
            ax.grid(True)
        draw_x, = axes[0].plot([],'-o')
        draw_y, = axes[1].plot([],'-o')
        draw_r, = axes[2].plot([],'-o')
        fig.canvas.draw()
        fig.show()

    # scan
    for ii, angle in enumerate(angles_1):
        angle_1 = angle
        angle_2 = angles_2[ii]
        if (angle_1 == start_angle):
            #move_axis_1(angle_1-move_back_1, axis=axis_1)
            #move_axis_2(angle_2-move_back_2, axis=axis_2)
            #move_axis_1(angle_1, axis=axis_1)
            #move_axis_2(angle_2, axis=axis_2)
            newport.corotate_axes(1, 2, angle_1-move_back_1, angle_2-move_back_2, axis_1=axis_1, axis_2=axis_2)
            #ctrl.corotate_axes(1, 2, angle_1, angle_2, axis_1=axis_1, axis_2=axis_2)
            time.sleep(2) # was set to 2 but I don't think this matters much
        newport.corotate_axes(1, 2, angle_1, angle_2, axis_1=axis_1, axis_2=axis_2)
        #move_axis_1(angle_1, axis=axis_1)
        #move_axis_2(angle_2, axis=axis_2)
        time.sleep(time_constant)

        # read lockin, rotators
        lockin_data = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index)
        angle_pos_1 = read_axis_1(axis=axis_1, print_flag=False)
        angle_pos_2 = read_axis_2(axis=axis_2, print_flag=False)

        # extract lockin data
        x = lockin_data['Demod x']
        y = lockin_data['Demod y']
        r = lockin_data['Demod r']

        # read additional motors
        measured_positions_dict = read_motors(mobj_measure_dict)
        measured_positions = [measured_positions_dict[m] for m in measure_motors]

        position = np.append(position, angle_pos_1)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # write to file
        if savefile:
            vars = [angle_pos_1, angle_pos_2]+list(lockin_data.values())+measured_positions
            append_data_to_file(fname, vars)

        # update plot
        if showplot == True:
            draw_x.set_data(position, demod_x)
            draw_y.set_data(position, demod_y)
            draw_r.set_data(position, demod_r)
            for ax in axes:
                ax.relim()
                ax.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

    # move motors back to original positions
    #move_axis_1(start_angle, axis=axis_1)
    #move_axis_2(start_angle+angle_offset, axis=axis_2)
    newport.corotate_axes(1, 2, start_angle, start_angle+angle_offset, axis_1=axis_1, axis_2=axis_2)

    # close motors
    if not passed_measure_motors:
        close_motors(mobj_measure_dict)

    return position, demod_x, demod_y, demod_r

def motor_scan(map_dict, filename_head=None, filename=None, showplot=True, time_constant=0.3, channel_index=1, R_channel_index=2, print_flag=False, savefile=True, metadata={}, daq_objs=None):
    '''
    utility to record lockin measurement as a function of motors specified by dictionary map_dict.
    '''

    # Lock-in Amplifier initialization
    if daq_objs is None:
        daq_objs = instrument_dict['zurich_lockin']['init']()
    read_lockin = instrument_dict['zurich_lockin']['read']
    lockin_header = list(read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index).keys())

    # capture motor information and initialize
    motors, mranges, mkwargs_dict = capture_motor_information(map_dict)
    measure_motors = list(set(meta_motors) - set(motors))
    mobj_dict = initialize_motors(motors)
    mobj_measure_dict = initialize_motors(measure_motors)

    # setup metadata - ie, for quick reference of starting state before measurement
    metadata = {**metadata, **generate_metadata({**mobj_dict, **mobj_measure_dict})}

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)
    #print(positions)

    # setup file for writing
    if savefile:
        append=False
        if filename==None:
            append=True
        filename_head, filename = generate_filename(filename_head, filename, 'motor_scan')
        if append==True:
            for m in motors:
                filename = filename+f'_{m}'
        fname = get_unique_filename(filename_head, filename)
        header_motors = [motor_dict[m]['name']+str(' measured') for m in measure_motors]
        header = [motor_dict[m]['name'] for m in motors]+[motor_dict[m]['name']+str(' measured') for m in motors]+lockin_header+header_motors
        write_file_header(fname, header, metadata)

    # move motors to start position, using move_back to handle initial case
    move_motors_to_start(mobj_dict, mkwargs_dict, positions, print_flag=print_flag)

    # setup measureables
    recorded_positions = [np.array([]) for i in range(len(map_dict))]
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # setup plots
    if showplot==True:
        if len(map_dict)==1:

            fig, axes = plt.subplots(3, 1, figsize=(8,10))
            x_label = motor_dict[motors[0]]['name']
            y_labels = ['Demod x', 'Demod y', 'R']
            for ii, ax in enumerate(axes):
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_labels[ii])
                ax.grid(True)
            draw_x, = axes[0].plot([],'-o')
            draw_y, = axes[1].plot([],'-o')
            draw_r, = axes[2].plot([],'-o')
            fig.canvas.draw()
            fig.show()

        elif len(map_dict)==2:

            # setup measureables
            xrange = mranges[0]
            yrange = mranges[1]
            x_num = len(xrange)
            y_num = len(yrange)
            demod_x0 = np.zeros((y_num, x_num))
            demod_y0 = np.zeros((y_num, x_num))
            demod_r0 = np.zeros((y_num, x_num))
            X_coor, Y_coor = np.meshgrid(xrange, yrange)
            extent=[xrange[0], xrange[-1], yrange[0], yrange[-1]]

            fig, axes = plt.subplots(3, 1, figsize=(8,10))
            x_label = motor_dict[motors[0]]['name']
            y_label = motor_dict[motors[1]]['name']
            for ii, ax in enumerate(axes):
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)
            mapx = axes[0].imshow(demod_x0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0, vmin=min(demod_x0.min(), -1e-12), vmax=max(demod_x0.max(), 1e-12)))
            mapy = axes[1].imshow(demod_y0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0, vmin=min(demod_y0.min(), -1e-12), vmax=max(demod_y0.max(), 1e-12)))
            mapr = axes[2].imshow(demod_r0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0, vmin=min(demod_r0.min(), -1e-12), vmax=max(demod_r0.max(), 1e-12)))
            fig.colorbar(mapx, ax=axes[0])
            fig.colorbar(mapy, ax=axes[1])
            fig.colorbar(mapr, ax=axes[2])
            fig.canvas.draw()
            fig.tight_layout()
            fig.show()
        else:
            print('Cannot plot scans that are greater than 2 dimensions.')

    # loop over positions, only moving a motor if its target position has changed.
    start_pos = positions[0]
    current_pos = start_pos
    num_pos = len(positions)
    for ii in tqdm(range(num_pos)):
        pos = positions[ii]

        # move motors to go back position if starting new raster - NEEDS WORK
        #if pos[0]!=current_pos[0] and pos[1]!=current_pos[1]:
        #    move_motors(mobj_dict, mkwargs_dict, current_pos, pos-10)

        # move motors if position has changed
        move_motors(mobj_dict, mkwargs_dict, current_pos, start_pos, pos, print_flag=print_flag)
        current_pos = pos

        # acquire data
        lockin_data = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index)

        # extract lockin data
        x = lockin_data['Demod x']
        y = lockin_data['Demod y']
        r = lockin_data['Demod r']

        # read actual motor positions
        measured_positions_dict = read_motors(mobj_dict)
        measured_positions = [measured_positions_dict[m] for m in motors]

        # read additional motors
        measured_motors_positions_dict = read_motors(mobj_measure_dict)
        measured_motors_positions = [measured_motors_positions_dict[m] for m in measure_motors]

        # update measurable
        for ii, p in enumerate(measured_positions):
            recorded_positions[ii] = np.append(recorded_positions[ii], p)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # add to file
        if savefile:
            append_data_to_file(fname, list(pos)+list(measured_positions)+list(lockin_data.values())+list(measured_motors_positions))

        # update plots
        if showplot==True:
            if len(map_dict)==1:
                draw_x.set_data(recorded_positions[0], demod_x)
                draw_y.set_data(recorded_positions[0], demod_y)
                draw_r.set_data(recorded_positions[0], demod_r)
                for ax in axes:
                    ax.relim()
                    ax.autoscale()
                fig.canvas.draw()
                fig.canvas.flush_events()

            elif len(map_dict)==2:
                time.sleep(0.01)
                length = len(demod_x)
                y_num0 = length//x_num
                x_num0 = length-y_num0*x_num
                demod_x0[:y_num0, :] = np.reshape(demod_x[:y_num0*x_num], (y_num0, x_num))
                demod_y0[:y_num0, :] = np.reshape(demod_y[:y_num0*x_num], (y_num0, x_num))
                demod_r0[:y_num0, :] = np.reshape(demod_r[:y_num0*x_num], (y_num0, x_num))
                if (y_num0 < y_num):
                    demod_x0[y_num0, :x_num0] = demod_x[y_num0*x_num:length]
                    demod_y0[y_num0, :x_num0] = demod_y[y_num0*x_num:length]
                    demod_r0[y_num0, :x_num0] = demod_r[y_num0*x_num:length]

                #print(f'x: {demod_x0.min()}, {demod_x0.max()}')
                #print(f'y: {demod_y0.min()}, {demod_y0.max()}')
                #print(f'r: {demod_r0.min()}, {demod_r0.max()}')
                #print(demod_x0)
                mapx.set_data(demod_x0)
                mapx.set_clim(vmin = min(demod_x0.min(),0), vmax = max(demod_x0.max(),0))
                mapy.set_data(demod_y0)
                mapy.set_clim(vmin = min(demod_y0.min(),0), vmax = max(demod_y0.max(),0))
                mapr.set_data(demod_r0)
                mapr.set_clim(vmin = min(demod_r0.min(),0), vmax = max(demod_r0.max(),0))
                fig.canvas.draw()
                fig.canvas.flush_events()

    # close motors
    close_motors(mobj_dict)
    close_motors(mobj_measure_dict)

def motor_scan_balance(map_dict, balance, balance_table=None, slope=0, tol=0, balance_channel=3, autobalance_flag=True, filename_head=None, filename=None, showplot=True, time_constant=0.3, channel_index=1, R_channel_index=4, print_flag=False, savefile=True, metadata={}, daq_objs=None):
    '''
    utility to record lockin measurement as a function of motors specified by dictionary map_dict.

    autobalance before each scan of last motor in map_dict according to balance table, which has the form of an xarray object with balance angle data over the various coordinates.

    only works if corotate_axes12 is in the map_dict
    '''

    # setup threading stuff
    moving_motors = LockedVar(False)

    # setup autobalancing thread
    def autobalance_thread(slope, tolerance, kwargs):
        #print('starting thread')
        time.sleep(3) # only trigger if move is long
        daq_objs = kwargs['daq_objs']
        zurich.set_zurich_acfilter(0,daq_objs=daq_objs)
        while moving_motors.locked_read()==True:
            autobalance_cont(slope, tolerance, **kwargs)
            #print('iterating thread')
            time.sleep(0.5)
        zurich.set_zurich_acfilter(0,daq_objs=daq_objs)
        #print('ending thread')

    def move_motors_thread(mobj_dict, mkwargs_dict, current_pos, start_pos, pos, kwargs):
        move_motors(mobj_dict, mkwargs_dict, current_pos, start_pos, pos, **kwargs)
        moving_motors.locked_update(False)

    bal_table_flag = True
    if balance_table is None:
        bal_table_flag=False

    # Lock-in Amplifier initialization
    if daq_objs is None:
        daq_objs = instrument_dict['zurich_lockin']['init']()
    read_lockin = instrument_dict['zurich_lockin']['read']
    lockin_header = list(read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index).keys())

    # capture motor information and initialize
    motors, mranges, mkwargs_dict = capture_motor_information(map_dict)
    measure_motors = list(set(meta_motors) - set(motors))
    mobj_dict = initialize_motors(motors)
    mobj_measure_dict = initialize_motors(measure_motors)

    # setup metadata - ie, for quick reference of starting state before measurement
    metadata = {**metadata, **generate_metadata({**mobj_dict, **mobj_measure_dict})}


    ''' must come up with a good test to check.
    if 'corotate_axes12' in motors:
        coords = motors.copy()
        coords.remove('corotate_axes12')
        coords = coords+['axis_1']
        coords = [motor_dict[c]['name'] for c in coords]
    else:
        coords = motors.copy()
        coords = [motor_dict[c]['name'] for c in coords]
    for c in list(balance_table.coords):
        if c not in coords:
            raise ValueError(f'balance table coordinates {list(balance_table.coords)} not in map_dict coordinate {coords}.')
    '''

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)
    #print(positions)

    # setup file for writing
    if savefile:
        append=False
        if filename==None:
            append=True
        filename_head, filename = generate_filename(filename_head, filename, 'motor_scan')
        if append==True:
            for m in motors:
                filename = filename+f'_{m}'
        fname = get_unique_filename(filename_head, filename)
        header_motors = [motor_dict[m]['name']+str(' measured') for m in measure_motors]
        header = [motor_dict[m]['name'] for m in motors]+[motor_dict[m]['name']+str(' measured') for m in motors]+lockin_header+header_motors+['Balance Angle (deg)']
        write_file_header(fname, header, metadata)

    # move motors to start position, using move_back to handle initial case
    move_motors_to_start(mobj_dict, mkwargs_dict, positions, print_flag=print_flag)

    # setup measureables
    recorded_positions = [np.array([]) for i in range(len(map_dict))]
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # setup plots
    if showplot==True:
        if len(map_dict)==1:

            fig, axes = plt.subplots(3, 1, figsize=(8,10))
            x_label = motor_dict[motors[0]]['name']
            y_labels = ['Demod x', 'Demod y', 'R']
            for ii, ax in enumerate(axes):
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_labels[ii])
                ax.grid(True)
            draw_x, = axes[0].plot([],'-o')
            draw_y, = axes[1].plot([],'-o')
            draw_r, = axes[2].plot([],'-o')
            fig.canvas.draw()
            fig.show()

        elif len(map_dict)==2:

            # setup measureables
            xrange = mranges[0]
            yrange = mranges[1]
            x_num = len(xrange)
            y_num = len(yrange)
            demod_x0 = np.zeros((y_num, x_num))
            demod_y0 = np.zeros((y_num, x_num))
            demod_r0 = np.zeros((y_num, x_num))
            X_coor, Y_coor = np.meshgrid(xrange, yrange)
            extent=[xrange[0], xrange[-1], yrange[0], yrange[-1]]

            fig, axes = plt.subplots(3, 1, figsize=(8,10))
            x_label = motor_dict[motors[0]]['name']
            y_label = motor_dict[motors[1]]['name']
            for ii, ax in enumerate(axes):
                ax.set_xlabel(x_label)
                ax.set_ylabel(y_label)
            mapx = axes[0].imshow(demod_x0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0, vmin=min(demod_x0.min(), -1e-12), vmax=max(demod_x0.max(), 1e-12)))
            mapy = axes[1].imshow(demod_y0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0, vmin=min(demod_y0.min(), -1e-12), vmax=max(demod_y0.max(), 1e-12)))
            mapr = axes[2].imshow(demod_r0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0, vmin=min(demod_r0.min(), -1e-12), vmax=max(demod_r0.max(), 1e-12)))
            fig.colorbar(mapx, ax=axes[0])
            fig.colorbar(mapy, ax=axes[1])
            fig.colorbar(mapr, ax=axes[2])
            fig.canvas.draw()
            fig.tight_layout()
            fig.show()
        else:
            print('Cannot plot scans that are greater than 2 dimensions.')

    # setup for autobalancing
    if balance!='all':
        balance_idx = motors.index(balance)
    else:
        balance_idx = 0
    if 'corotate_axes12' in motors:
        axis_1, axis_2 = mobj_dict['corotate_axes12']
    else:
        if 'axis_1' in motors:
            axis_1 = mobj_dict['axis_1']
        else:
            axis_1 = motor_dict['axis_1']['init']()
        if 'axis_2' in motors:
            axis_2 = mobj_dict['axis_2']
        else:
            axis_2 = motor_dict['axis_2']['init']()
    move_axis_1 = motor_dict['axis_1']['move']
    move_axis_2 = motor_dict['axis_2']['move']
    read_axis_1 = motor_dict['axis_1']['read']
    read_axis_2 = motor_dict['axis_2']['read']
    #corotate_axes_idx = motors.index('corotate_axes12')
    autobalkwrag_dict = {'daq_objs':daq_objs, 'axis_1':axis_1, 'axis_2':axis_2, 'channel_index':balance_channel}

    # loop over positions, only moving a motor if its target position has changed.
    start_pos = positions[0]
    current_pos = start_pos
    num_pos = len(positions)
    if 'corotate_axes12' in motors:
        bal_angle=mkwargs_dict['corotate_axes12']['bal_angle']
    else:
        bal_angle = read_axis_2()
    for ii in tqdm(range(num_pos)):
        pos = positions[ii]

        # move motors to go back position if starting new raster - NEEDS WORK
        #if pos[0]!=current_pos[0] and pos[1]!=current_pos[1]:
        #    move_motors(mobj_dict, mkwargs_dict, current_pos, pos-10)

        # move motors if position has changed
        bal_thread = threading.Thread(target=autobalance_thread, args=(slope, tol, autobalkwrag_dict))
        move_thread = threading.Thread(target=move_motors_thread, args=(mobj_dict, mkwargs_dict, current_pos, start_pos, pos, {'print_flag':print_flag}))
        moving_motors.locked_update(True)
        #bal_thread.start()
        move_thread.start()
        current_pos = pos
        move_thread.join()
        #bal_thread.join()

        if (current_pos[balance_idx] == start_pos[balance_idx]) or balance=='corotate_axes12' or balance=='all':
            if bal_table_flag==True:
                pos_dict = {}
                for ii, p in enumerate(pos):
                    m = motors[ii]
                    if m == 'corotate_axes12':
                        m = 'axis_1'
                        axis_1_pos = p
                    name = motor_dict[m]['name']
                    pos_dict[name] = p
                if 'axis_1' not in list(pos_dict.keys()):
                    name = motor_dict['axis_1']['name']
                    pos_dict[name] = read_axis_1()
                if list(balance_table.coords) == ['Corotation Axes (deg)']: # handle bal table with just angles
                    pos_dict = {'Corotation Axes (deg)':pos_dict['Angle 1 (deg)']}
                    bal_angle_approx, slope, tol = interp_balance_angle(pos_dict, balance_table)
                else:
                    bal_angle_approx, slope, tol = interp_balance_angle(pos_dict, balance_table)
                move_axis_2(axis_1_pos+bal_angle_approx, axis_2)
            if autobalance_flag==True:
                zurich.set_zurich_acfilter(0,daq_objs=daq_objs)
                bal_angle = autobalance(slope, tol, **autobalkwrag_dict)
                zurich.set_zurich_acfilter(1,daq_objs=daq_objs)
            elif bal_table_flag==True:
                bal_angle = bal_angle_approx
            if 'corotate_axes12' in motors:
                mkwargs_dict['corotate_axes12']['bal_angle'] = bal_angle

        # acquire data
        lockin_data = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index)

        # extract lockin data
        x = lockin_data['Demod x']
        y = lockin_data['Demod y']
        r = lockin_data['Demod r']

        # read actual motor positions
        measured_positions_dict = read_motors(mobj_dict)
        measured_positions = [measured_positions_dict[m] for m in motors]

        # read additional motors
        measured_motors_positions_dict = read_motors(mobj_measure_dict)
        measured_motors_positions = [measured_motors_positions_dict[m] for m in measure_motors]

        # update measurable
        for ii, p in enumerate(measured_positions):
            recorded_positions[ii] = np.append(recorded_positions[ii], p)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # add to file
        if savefile:
            append_data_to_file(fname, list(pos)+list(measured_positions)+list(lockin_data.values())+list(measured_motors_positions)+[bal_angle])

        # update plots
        if showplot==True:
            if len(map_dict)==1:
                draw_x.set_data(recorded_positions[0], demod_x)
                draw_y.set_data(recorded_positions[0], demod_y)
                draw_r.set_data(recorded_positions[0], demod_r)
                for ax in axes:
                    ax.relim()
                    ax.autoscale()
                fig.canvas.draw()
                fig.canvas.flush_events()

            elif len(map_dict)==2:
                time.sleep(0.01)
                length = len(demod_x)
                y_num0 = length//x_num
                x_num0 = length-y_num0*x_num
                demod_x0[:y_num0, :] = np.reshape(demod_x[:y_num0*x_num], (y_num0, x_num))
                demod_y0[:y_num0, :] = np.reshape(demod_y[:y_num0*x_num], (y_num0, x_num))
                demod_r0[:y_num0, :] = np.reshape(demod_r[:y_num0*x_num], (y_num0, x_num))
                if (y_num0 < y_num):
                    demod_x0[y_num0, :x_num0] = demod_x[y_num0*x_num:length]
                    demod_y0[y_num0, :x_num0] = demod_y[y_num0*x_num:length]
                    demod_r0[y_num0, :x_num0] = demod_r[y_num0*x_num:length]

                #print(f'x: {demod_x0.min()}, {demod_x0.max()}')
                #print(f'y: {demod_y0.min()}, {demod_y0.max()}')
                #print(f'r: {demod_r0.min()}, {demod_r0.max()}')
                #print(demod_x0)
                mapx.set_data(demod_x0)
                mapx.set_clim(vmin = min(demod_x0.min(),0), vmax = max(demod_x0.max(),0))
                mapy.set_data(demod_y0)
                mapy.set_clim(vmin = min(demod_y0.min(),0), vmax = max(demod_y0.max(),0))
                mapr.set_data(demod_r0)
                mapr.set_clim(vmin = min(demod_r0.min(),0), vmax = max(demod_r0.max(),0))
                fig.canvas.draw()
                fig.canvas.flush_events()

    # close motors
    close_motors(mobj_dict)
    close_motors(mobj_measure_dict)

def rotate_map(map_dict, start_angle, end_angle, step_size, filename_head=None, filename=None, axis_index=1, showplot=False, time_constant=0.3, channel_index=1, R_channel_index=4, daq_objs=None, print_flag=False, savefile=True, metadata={}):

    # Lock-in Amplifier initialization
    daq_objs = instrument_dict['zurich_lockin']['init']()

    # initialize rotation axes
    axis_1 = motor_dict['axis_1']['init']()
    axis_2 = motor_dict['axis_2']['init']()

    # capture motor information and initialize
    motors, mranges, mkwargs_dict = capture_motor_information(map_dict)
    measure_motors = list(set(meta_motors) - set(motors) - set(['axis_1', 'axis_2']))
    mobj_dict = initialize_motors(motors)
    mobj_measure_dict = initialize_motors(measure_motors)

    # setup metadata - ie, for quick reference of starting state before measurement
    metadata = {**metadata, **generate_metadata({**mobj_dict, **mobj_measure_dict, **{'axis_1':axis_1, 'axis_2':axis_2}})}

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)

    # move motors to start position, using move_back to handle initial case
    move_motors_to_start(mobj_dict, mkwargs_dict, positions, print_flag=print_flag)

    # setup filenames
    filename_head, filename = generate_filename(filename_head, filename, 'rotate_map')

    # loop over positions, only moving a motor if its target position has changed.
    start_pos = positions[0]
    current_pos = start_pos
    num_pos = len(positions)
    for ii in tqdm(range(num_pos)):
        pos = positions[ii]

        # move motors if position has changed
        move_motors(mobj_dict, mkwargs_dict, current_pos, start_pos, pos, print_flag=print_flag)

        # setup each filename
        expanded_filename = filename
        for ii, m in enumerate(motors):
            p = pos[ii]
            expanded_filename = expanded_filename+f'_{m}{p}'

        # scan
        rotate_scan(start_angle, end_angle, step_size, filename_head=filename_head, filename=expanded_filename, axis_index=axis_index, mobj_measure_dict={**mobj_dict, **mobj_measure_dict}, metadata=metadata, showplot=showplot, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2, savefile=savefile, override_metadata=True)

        current_pos = pos

    # close motors
    close_motors(mobj_dict)
    close_motors(mobj_measure_dict)

def corotate_map(map_dict, start_angle, end_angle, step_size, angle_offset, rate_axis_2=1, filename_head=None, filename=None, measure_motors=[], showplot=False, time_constant=0.3, channel_index=1, R_channel_index=4, print_flag=False, savefile=True, metadata={}):
    '''
    Takes a corotation scan at each point in a map specified by dictionary map_dict, which entries of the form 'axis':(start, end, step_size, kwargs), where kwargs is a dictionary of key/value pairs appropriate for each motor 'move' function. For example, a temperature map might take the following map dictionary:

    map_dict = {'temp':(10,20,1,{'tolerance':0.01, 'wait_time':30})}

    '''

    # Lock-in Amplifier initialization
    daq_objs = instrument_dict['zurich_lockin']['init']()

    # initialize rotation axes
    axis_1 = motor_dict['axis_1']['init']()
    axis_2 = motor_dict['axis_2']['init']()

    # capture motor information and initialize
    motors, mranges, mkwargs_dict = capture_motor_information(map_dict)
    measure_motors = list(set(meta_motors) - set(motors) - set(['axis_1', 'axis_2']))
    mobj_dict = initialize_motors(motors)
    mobj_measure_dict = initialize_motors(measure_motors)

    # setup metadata - ie, for quick reference of starting state before measurement
    metadata = {**metadata, **generate_metadata({**mobj_dict, **mobj_measure_dict, **{'axis_1':axis_1, 'axis_2':axis_2}})}

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)

    # move motors to start position, using move_back to handle initial case
    move_motors_to_start(mobj_dict, mkwargs_dict, positions, print_flag=print_flag)

    # setup filenames
    filename_head, filename = generate_filename(filename_head, filename, 'corotate_map')

    # loop over positions, only moving a motor if its target position has changed.
    start_pos = positions[0]
    current_pos = start_pos
    num_pos = len(positions)
    for ii in tqdm(range(num_pos)):
        pos = positions[ii]

        # move motors if position has changed, accounting for move_back
        move_motors(mobj_dict, mkwargs_dict, current_pos, start_pos, pos, print_flag=print_flag)

        # setup each filename
        expanded_filename = filename
        for ii, m in enumerate(motors):
            p = pos[ii]
            if m == 'temp':
                p = round(p,5)
            expanded_filename = expanded_filename+f'_{m}{p}'

        # scan
        corotate_scan(start_angle, end_angle, step_size, angle_offset, rate_axis_2=rate_axis_2, filename_head=filename_head, filename=expanded_filename, mobj_measure_dict={**mobj_dict, **mobj_measure_dict}, metadata=metadata, showplot=showplot, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2, savefile=savefile, override_metadata=True)

        current_pos = pos

    # close motors
    close_motors(mobj_dict)
    close_motors(mobj_measure_dict)

def monitor_motor(motor, filename_head=None, filename=None, savefile=False, showplot=True):
    '''
    Utility to measure the value of a motor as function of time.
    '''
    # setup file for writing
    if savefile:
        append=False
        if filename==None:
            append=True
        filename_head, filename = generate_filename(filename_head, filename, 'monitor_motor')
        if append==True:
                filename = filename+f'_{motor}'
        fname = get_unique_filename(filename_head, filename)
        header = ['Time (s)', motor_dict[motor]['name']]
        write_file_header(fname, header)

    # setup locked variable
    run = LockedVar(True)

    # setup thread function to stop function via user input
    def get_user_input(run_var):
        answer = input('Press any key to stop:')
        run_var.locked_update(False)

    mobj_dict = initialize_motors([motor])

    time_vect = []
    motor_vals = []

    if showplot==True:
        fig, ax = plt.subplots(1, 1, figsize=(6,4))
        x_label = 'Time (s)'
        y_label = motor_dict[motor]['name']
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.grid(True)
        draw_motor, = ax.plot([],'-o')
        fig.canvas.draw()
        fig.show()

    # start user input thread
    user_input_thread = threading.Thread(target=get_user_input, args=(run,))
    user_input_thread.start()

    t0 = time.time()
    while run.locked_read():

        measured_position = read_motors(mobj_dict)[motor]
        t = time.time() - t0

        time_vect.append(t)
        motor_vals.append(measured_position)

        # add to file
        if savefile:
            append_data_to_file(fname, [t, measured_position])

        # update plots
        if showplot==True:
            draw_motor.set_data(time_vect, motor_vals)
            ax.relim()
            ax.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()
        time.sleep(0.1)

    user_input_thread.join()

    # close motors
    close_motors(mobj_dict)

#########################
### Balancing Methods ###
#########################

def find_balance_angle_macro(start_angle, end_angle, step_size, step_size_fine=0.2, window=3, axis_index=2, balance_at=0, channel_index=1, time_constant=0.3, daq_objs=None, axis_1=None, axis_2=None):

    # initialize zurich lockin and setup read function
    if daq_objs==None:
        init_func = instrument_dict['zurich_lockin']['init']
        daq_objs = init_func()
    read_lockin = instrument_dict['zurich_lockin']['read']

    # initialize axes and setup move functions
    if axis_1 == None:
        init_func = motor_dict['axis_1']['init']
        axis_1 = init_func()
    if axis_2 == None:
        init_func = motor_dict['axis_2']['init']
        axis_2 = init_func()

    # initialize axes funcs
    move_axis_1 = motor_dict['axis_1']['move']
    move_axis_2 = motor_dict['axis_2']['move']
    move_back_1 = motor_dict['axis_1']['move_back']
    move_back_2 = motor_dict['axis_2']['move_back']

    # move both motors to 0
    move_axis_1(balance_at-move_back_1)
    move_axis_2(balance_at-move_back_2)
    move_axis_1(balance_at)
    move_axis_2(balance_at)

    coarse_scan = rotate_scan(start_angle+balance_at, end_angle+balance_at, step_size, channel_index=channel_index, time_constant=time_constant, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2, override_metadata=True, savefile=False, axis_index=2)

    fitf = lambda x, a, phi, c: a*np.cos(4*(2*np.pi/360)*(x-phi)) + c
    popt, pcov = opt.curve_fit(fitf, coarse_scan[0], coarse_scan[1], p0=[np.max(coarse_scan[1]), balance_at, 0], bounds=([0,-np.inf,-np.inf], [np.inf,np.inf,np.inf]))

    fig, ax = plt.subplots(1)
    xvect = np.linspace(start_angle+balance_at,end_angle+balance_at,1000)
    ys = [fitf(i, *popt) for i in xvect]
    ax.plot(coarse_scan[0], coarse_scan[1], 'o')
    ax.plot(xvect, ys, '-', color='black')

    a_opt = popt[0]
    phi_opt = popt[1]
    c_opt = popt[2]
    bal_angle_approx = (180/np.pi)*np.arccos(-c_opt/a_opt)/4 + phi_opt - balance_at
    print(bal_angle_approx)
    offset=c_opt

    bal_angle, slope, tol = find_balance_angle(bal_angle_approx-window, bal_angle_approx+window, step_size_fine, balance_at=balance_at, offset=offset, go_to_balance_angle=True, axis_index=2, channel_index=channel_index, time_constant=time_constant, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2)

    return bal_angle, slope, tol

def find_balance_angle(start_angle, end_angle, step_size, balance_at=0, offset=0, go_to_balance_angle=True, axis_index=2, channel_index=1, time_constant=0.3, R_channel_index=2, daq_objs=None, axis_1=None, axis_2=None):
    '''
    Assuming we are measuring in DC mode above a transition or on GaAs, carries out a rotate_scan. Find angle by carrying out a linear fit, such that the angle range should be taken to be very small.

    By default, moves axis 2 to find balance angle.

    If go_to_balance_angle is set to true, moves stages to
    '''

    # initialize zurich lockin and setup read function
    if daq_objs==None:
        init_func = instrument_dict['zurich_lockin']['init']
        daq_objs = init_func()
    read_lockin = instrument_dict['zurich_lockin']['read']

    # initialize axes and setup move functions
    if axis_1 == None:
        init_func = motor_dict['axis_1']['init']
        axis_1 = init_func()
    if axis_2 == None:
        init_func = motor_dict['axis_2']['init']
        axis_2 = init_func()

    # initialize axes funcs
    move_axis_1 = motor_dict['axis_1']['move']
    move_axis_2 = motor_dict['axis_2']['move']
    move_back_1 = motor_dict['axis_1']['move_back']
    move_back_2 = motor_dict['axis_2']['move_back']

    # move both motors to 0
    move_axis_1(balance_at-move_back_1)
    move_axis_2(balance_at-move_back_2)
    move_axis_1(balance_at)
    move_axis_2(balance_at)

    positions, demod_x, demod_y, demod_r = rotate_scan(balance_at+start_angle, balance_at+end_angle, step_size, axis_index=axis_index, channel_index=channel_index, time_constant=time_constant, R_channel_index=R_channel_index, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2, override_metadata=True, savefile=False)

    # linear fit
    fitf = lambda x, m, b: m*x+b
    popt, pcov = opt.curve_fit(fitf, positions, demod_x-offset)
    balance_angle = -popt[1]/popt[0]
    slope = popt[0]
    tol = np.abs(np.min(demod_x)*1.5)
    angles_vect = np.linspace(start_angle+balance_at, balance_at+end_angle, 1000)
    fit = popt[0]*angles_vect + popt[1]

    # display result
    fig, ax = plt.subplots(1)
    ax.set_ylabel('Angle (deg)')
    ax.set_ylabel('Demod X')
    ax.plot(positions, demod_x-offset, 'o', ms=5, color='blue')
    ax.plot(angles_vect, fit, '-', color='black')

    # set balance angle relative to 0
    balance_angle = balance_angle - balance_at

    if go_to_balance_angle: # is there a good way to automatically re-zero zero?
        if axis_index==1:
            move_axis_1(balance_angle+balance_at)
        else:
            move_axis_2(balance_angle+balance_at)

    print(f'Balance angle: {balance_angle}')
    return balance_angle, slope, tol

def autobalance_cont(slope, tolerance=None, daq_objs=None, axis_1=None, axis_2=None, balance_at=None, offset=0, channel_index=1, time_constant=0.3, print_flag=True, lockin_index=0):

    # initialize lockin
    if daq_objs==None:
        daq_objs = instrument_dict['zurich_lockin']['init']()
    else:
        pass
    read_lockin = instrument_dict['zurich_lockin']['read']

    # initialize axes and setup axis move functions
    if axis_1 == None:
        init_func = motor_dict['axis_1']['init']
        axis_1 = init_func()
    if axis_2 == None:
        init_func = motor_dict['axis_2']['init']
        axis_2 = init_func()
    move_axis_1 = motor_dict['axis_1']['move']
    move_axis_2 = motor_dict['axis_2']['move']
    read_axis_1 = motor_dict['axis_1']['read']
    read_axis_2 = motor_dict['axis_2']['read']
    move_back_1 = motor_dict['axis_1']['move_back']
    move_back_2 = motor_dict['axis_2']['move_back']

    # move axis 1 to balance_at position if value given
    if balance_at!=None:
        move_axis_1(balance_at-move_back_1, axis=axis_1)
        move_axis_1(balance_at, axis=axis_1)
    else:
        balance_at = read_axis_1(axis=axis_1)

    curr_pos = read_axis_2(axis=axis_2)
    pid_signal = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index)[lockin_index]
    new_pos = curr_pos-(1/slope)*(pid_signal-offset)
    move_axis_2(new_pos, axis=axis_2)

def autobalance(slope, tolerance, daq_objs=None, axis_1=None, axis_2=None, balance_at=None, offset=0, channel_index=1, time_constant=0.3, print_flag=True):

    '''
    balance lockin at specified balance_at angle on the fly. Axis 1 is taken to be held fixed at balance_at and axis 2 is moved to balance PID (photodiode).

    in the future this could be designed more robustly with a PID control loop
    '''

    # initialize lockin
    if daq_objs==None:
        daq_objs = instrument_dict['zurich_lockin']['init']()
    else:
        pass
    read_lockin = instrument_dict['zurich_lockin']['read']

    # initialize axes and setup axis move functions
    if axis_1 == None:
        init_func = motor_dict['axis_1']['init']
        axis_1 = init_func()
    if axis_2 == None:
        init_func = motor_dict['axis_2']['init']
        axis_2 = init_func()
    move_axis_1 = motor_dict['axis_1']['move']
    move_axis_2 = motor_dict['axis_2']['move']
    read_axis_1 = motor_dict['axis_1']['read']
    read_axis_2 = motor_dict['axis_2']['read']
    move_back_1 = motor_dict['axis_1']['move_back']
    move_back_2 = motor_dict['axis_2']['move_back']

    # move axis 1 to balance_at position if value given
    if balance_at!=None:
        move_axis_1(balance_at-move_back_1, axis=axis_1)
        move_axis_1(balance_at, axis=axis_1)
    else:
        balance_at = read_axis_1(axis=axis_1)

    pid_signal = 10000
    curr_pos = read_axis_2(axis=axis_2)
    while (np.abs(pid_signal-offset)>tolerance):
        #print(curr_pos)
        time.sleep(time_constant*4)
        pid_signal = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index)[f'Demod {channel_index} x']
        new_pos = (curr_pos-(1/slope)*(pid_signal-offset))%360
        move_axis_2(new_pos, axis=axis_2)
        curr_pos = new_pos
    if print_flag:
        print(f'Balanced PID at {balance_at}. Balance angle: {curr_pos-balance_at}.')
    return curr_pos-balance_at

def create_balance_table(map_dict, start_angle=0, end_angle=90, step_size=5, step_size_fine=0.2, window=3, offset=0, filename_head=None, filename=None, channel_index=1, time_constant=0.3, daq_objs=None, axis_1=None, axis_2=None, print_flag=False, save_metadata=True):
    '''
    create balance table for axis_1 and any other parameters in map_dict.
    '''

    # setup metadata
    metadata = generate_metadata()

    # capture motor information and initialize
    motors, mranges, mkwargs_dict = capture_motor_information(map_dict)
    mobj_dict = initialize_motors(motors)
    #if 'axis_1' not in motors:
    #    raise ValueError('axis_1 must be in the map_dict, even if making a balance table at a single angle')
    #if 'axis_2' in motors:
    #    raise ValueError('axis 2 cannot be in the map_dict')

    # initialize rotation axes
    if 'axis_1' in motors:
        axis_1 = mobj_dict['axis_1']
    else:
        axis_1 = motor_dict['axis_1']['init']()
    if 'axis_2' in motors:
        axis_2 = mobj_dict['axis_2']
    else:
        axis_2 = motor_dict['axis_2']['init']()

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)
    #print(positions)

    # setup file for writing
    filename_head, filename = generate_filename(filename_head, filename, 'balance_table')
    for m in motors:
        filename = filename+f'_{m}'
    fname = get_unique_filename(filename_head, filename)
    header = [motor_dict[m]['name'] for m in motors]+['Balance Angle (deg)', 'Slope', 'Tolerance']
    write_file_header(fname, header, metadata)

    # move motors to start position, using move_back to handle initial case
    move_motors_to_start(mobj_dict, mkwargs_dict, positions, print_flag=print_flag)

    # loop over positions, only moving a motor if its target position has changed.
    if 'axis_1' in motors:
        angle1_idx = motors.index('axis_1')
    else:
        balance_at = motor_dict['axis_1']['read']()
    start_pos = positions[0]
    current_pos = start_pos
    num_pos = len(positions)
    for ii in tqdm(range(num_pos)):
        pos = positions[ii]

        # move motors to go back position if starting new raster - NEEDS WORK
        #if pos[0]!=current_pos[0] and pos[1]!=current_pos[1]:
        #    move_motors(mobj_dict, mkwargs_dict, current_pos, pos-10)

        # move motors if position has changed
        move_motors(mobj_dict, mkwargs_dict, current_pos, start_pos, pos, print_flag=print_flag)
        current_pos = pos

        # find balance angle
        if 'axis_1' in motors:
            balance_at = pos[angle1_idx]
        bal_angle, slope, tol = find_balance_angle_macro(start_angle, end_angle, step_size, step_size_fine=step_size_fine, window=window, axis_index=2, balance_at=balance_at, offset=offset, channel_index=channel_index, time_constant=time_constant, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2)

        # add to file
        append_data_to_file(fname, list(pos)+[bal_angle, slope, tol])

    # close motors
    close_motors(mobj_dict)

    # process the results using orensteinanalysis and pickle
    motor_names = [motor_dict[m]['name'] for m in motors]
    obj = loader.load_measurement(fname)
    obj  = process.reshape(obj, motor_names)
    obj_filename = f'{fname[:-4]}_obj.cdf'
    obj.to_netcdf(obj_filename)


    return obj

def create_balance_table_quick(slope, tol, start_angle, end_angle, step_size, bal_angle_init, offset=0, filename_head=None, filename=None, channel_index=1, time_constant=0.3, daq_objs=None, axis_1=None, axis_2=None, print_flag=False, save_metadata=True):
    '''
    create balance table for axis_1 and any other parameters in map_dict.
    '''

    # setup metadata
    metadata = generate_metadata()

    daq_objs = instrument_dict['zurich_lockin']['init']()

    # initialize rotation axes
    axis_1 = motor_dict['axis_1']['init']()
    axis_2 = motor_dict['axis_2']['init']()

    # setup file for writing
    filename_head, filename = generate_filename(filename_head, filename, 'balance_table')
    fname = get_unique_filename(filename_head, filename)
    header = ['Corotation Axes (deg)']+['Balance Angle (deg)', 'Slope', 'Tolerance']
    write_file_header(fname, header, metadata)

    positions = get_motor_range(start_angle, end_angle, step_size)

    # move motors to start position
    motor_dict['axis_1']['move'](0, axis_1)
    motor_dict['axis_2']['move'](0+bal_angle_init, axis_2)

    bal_angle = bal_angle_init
    for ii in tqdm(range(len(positions))):
        pos = positions[ii]

        # move motors
        motor_dict['axis_1']['move'](pos, axis_1)
        motor_dict['axis_2']['move'](pos+bal_angle, axis_2)

        # find balance angle
        bal_angle = autobalance(slope, tol, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2, channel_index=channel_index)
        # add to file
        append_data_to_file(fname,[pos]+[bal_angle, slope, tol])

    # process the results using orensteinanalysis and pickle
    obj = loader.load_measurement(fname)
    obj  = process.reshape(obj, ['Corotation Axes (deg)'])
    obj_filename = f'{fname[:-4]}_obj.cdf'
    obj.to_netcdf(obj_filename)

    return obj

def interp_balance_angle(pos_dict, balance_table):
    '''
    given a balance table and a dictionary giving positions within that balance table, interpolate the value of balance angle, slope, and tolerance.
    '''

    bal_angle_data = balance_table['Balance Angle (deg)'].data
    slope_data = balance_table['Slope'].data
    tol_data = balance_table['Tolerance'].data
    motors = list(balance_table.dims)
    points = tuple([balance_table[m].data for m in motors])

    pos = np.array([pos_dict[m] for m in motors])
    bal_angle_interp = interp.interpn(points, bal_angle_data, pos)
    slope_interp = interp.interpn(points, slope_data, pos)
    tol_interp = interp.interpn(points, tol_data, pos)

    return bal_angle_interp[0], slope_interp[0], tol_interp[0]

###################
### Lab Helpers ###
###################

def align_delay_stage(wait_time=5, range=(-125,125)):

    # initialize delay stage:
    delay_stage = motor_dict['delay_stage']['init']()
    move = motor_dict['delay_stage']['move']

    # setup locked variable
    run = True
    lock = threading.Lock()

    # setup thread func
    def osc_delay_stage(wait):
        run_cond = True
        n=0
        while run_cond:
            if n!=0:
                time.sleep(wait)
            if n==0:
                n=1
            move(range[0], delay_stage)
            time.sleep(wait)
            move(range[1], delay_stage)
            with lock: # read locked variable
                run_cond = run
                # print(run_cond)

    # start thread
    move_ds_thread = threading.Thread(target=osc_delay_stage, args=(wait_time,))
    move_ds_thread.start()

    # ask to stop
    answer = input('Press any key to stop:')

    # shutdown
    with lock:
        run = False
    move_ds_thread.join()

def list_motors():
    '''
    walk through motors in motor_dict and return those that can initialize properly
    '''

    m_list = []
    for m in motor_dict.keys():
        try:
            obj = motor_dict[m]['init']()
            m_list.append(m)
            motor_dict[m]['close'](obj)
        except:
            print(f'motor {m} not able to initialize')

    return m_list

###########################
### Strain Cell Methods ### # perhaps move these to another file?
###########################

def measure_strain_cell_capacitor(fname, sc, num_points=1000, dt_min=0.1):
    '''
    mesaures the
    '''
    t = np.zeros(num_points)
    cap = np.zeros(num_points)
    t0 = time.time()
    i = 0
    t_old = t0
    while i < num_points:
        t_new = time.time()
        if t_new - t_old > dt_min:
            c = sc.get_cap()
            cap[i] = c
            t[i] = t_new - t0
            i = i + 1
            t_old = t_new
    save_data_to_file(fname, np.transpose([t, cap]), ['Time', 'Cap'])

def strain_cell_temperature_calibration(fname1, fname2, filename_head, sc, cryo, temps, lakeshore_stability, cap_stability, cryo_stability_high, cryo_stability_low, mode, wait_time=1):
    '''
    runs a cooldown and warmup of Montana CryoAdvance and monitors platform temperature vs strain cell capacitance. The strain cell should be loaded with a titanium dummy sample. For consistency, block temperature should be measured NOT sample temperature, which will change from setup to setup.

    args:
        - fname1:           file path to a file to log continuosly during a ramp
        - fname2:           file path to a file to log only at times when both temperature and capacitor are stable
        - sc:               StrainClient object
        - cryo:             Montana CryoCore object
        - temps:            list of temperatures. Note that Montana seems to like integer values.
        - mode:             1 for keeping cryostat at base and using heater and 2 for changing cold plate temp
    '''
    BASE_TEMP=12.0 # needs to actually reach base to read as stable in Montana software, so make a bit higher than base.
    filename1 = filename_head+'\\'+fname1+'.dat'
    filename2 = filename_head+'\\'+fname2+'.dat'
    with open(filename1, 'a') as f1:
        f1.write('Time' + '\t' + 'Setpoint Temperature (K)' + '\t' + 'Platform Temperature (K)' + '\t' + 'Lakeshore Temperature (K)' + '\t' + 'Capacitance' + '\n')
    with open(filename2, 'a') as f2:
        f2.write('Platform Temperature (K)' + '\t' + 'Lakeshore Temperature (K)' + '\t' + 'Capacitance' + '\n')

    t0 = time.time()
    cryo.cooldown()
    print('Cooling down cryostat.')

    if mode==1: # cooldown to base temperature and then set temperature with lakeshore, waiting for both lakeshore and capacitance to stabilize before moving on.
        print('Measuring in mode 1.')
        setpoints = np.sort(temps) # measure from base to high temp.
        target_stability = cryo_stability_low
        cryo.set_platform_target_temperature(BASE_TEMP)
        cryo.set_platform_stability_target(target_stability)
        ctrl.set_temperature(0)
        ctrl.set_lakeshore_range(0)
        while True:
            time.sleep(wait_time)
            stability_ok, is_stable = cryo.get_platform_temperature_stable()
            if is_stable:
                print(f'Stabilized platform to base temperature.')
                break
        cryo.set_platform_target_temperature(0)
        for sp in setpoints:
            if sp < 9:
                ctrl.set_lakeshore_range(1)
            elif sp < 11 and sp > 9:
                ctrl.set_lakeshore_range(2)
            else:
                ctrl.set_lakeshore_range(3)
            print(f'Setting setpoint to {sp} K')
            ctrl.set_temperature(sp)
            lakeshore_temps = []
            while True:
                time.sleep(wait_time)
                lakeshore_temps.append(ctrl.read_temperature())
                with open(filename1, 'a') as f1:
                    f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                         str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                         str(format(float(lakeshore_temps[-1]), '.5f')) + '\t' +
                         str(format(float(sc.get_cap()), '.5f')) + '\n')
                if len(lakeshore_temps) > 120:
                    mean = np.mean(np.asarray(lakeshore_temps[-120:]))
                    std = np.std(np.asarray(lakeshore_temps[-120:]))
                    if ((std < lakeshore_stability) or (len(lakeshore_temps) > 7200)) and ((mean > sp-0.01) and (mean < sp+0.01)):
                        if std < lakeshore_stability:
                            print(f'Stabilized Lakeshore temperature at {mean} K')
                        if len(lakeshore_temps) > 7200:
                            print('Exceeded maximum soak time')
                        print(f'Lakeshore noise: {std}')
                        caps = []
                        while True:
                            time.sleep(wait_time)
                            caps.append(sc.get_cap())
                            lakeshore_temps.append(ctrl.read_temperature())
                            with open(filename1, 'a') as f1:
                                f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                                     str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                     str(format(float(lakeshore_temps[-1]), '.5f')) + '\t' +
                                     str(format(float(caps[-1]), '.5f')) + '\n')
                            if len(caps)>300: # mininum soak time of 5 minutes
                                if (np.std(np.asarray(caps[-300:])) < cap_stability) or (len(lakeshore_temps) > 7200 and len(caps) > 600):
                                    print('STD: '+str(np.std(np.asarray(caps[-300:]))))
                                    if (np.std(np.asarray(caps[-300:])) < cap_stability):
                                        print('Stdev below accepted value')
                                    elif len(caps) > 7200: # maximum soak time of 2 hours
                                        print('Exceeded maximum soak time')
                                    print(f'Stabilized capacitance measurement, writing to file')
                                    with open(filename2, 'a') as f2:
                                        f2.write(str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                             str(format(np.mean(lakeshore_temps[-15:]), '.5f')) + '\t' +
                                             str(format(np.mean(caps[-15:]), '.5f')) + '\n')
                                    break
                                elif (np.std(np.asarray(caps[-300:])) < cap_stability*1.1) and (len(caps) > 3600):
                                    print('STD: '+str(np.std(np.asarray(caps[-300:]))))
                                    print('Stdev within 10 percent of accepted value after 1 hour')
                                    print(f'Stabilized capacitance measurement, writing to file')
                                    with open(filename2, 'a') as f2:
                                        f2.write(str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                             str(format(float(ctrl.read_temperature()), '.5f')) + '\t' +
                                             str(format(np.mean(caps[-15:]), '.5f')) + '\n')
                                    break
                        break

    if mode==2: # on cooldown, set cryostat target temperature and wait for lakeshore and capacitance - this still needs some work to work around stability issues at 10-12K by switching the stability criteria for both cryostat and lakeshore.
        print('Measuring in mode 1.')
        setpoint = np.flip(np.sort(temps)) # measure from high temp to low temp.
        for sp in setpoints:
            print(f'Setting setpoint to {sp} K')
            if sp >= 10:
                target_stability = cryo_stability_high
            else:
                target_stability = cryo_stability_low
            cryo.set_platform_target_temperature(int(sp))
            cryo.set_platform_stability_target(target_stability)
            while True:
                time.sleep(wait_time)
                with open(filename1, 'a') as f1:
                    f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                         str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                         str(format(float(ctrl.read_temperature()), '.5f')) + '\t' + str(format(float(sc.get_cap()), '.5f')) + '\n')
                stability_ok, is_stable = cryo.get_platform_temperature_stable()
                if is_stable:
                    print(f'Stabilized platform temperature at {cryo.get_platform_temperature()[1]} K')
                    lakeshore_temps = []
                    while True:
                        time.sleep(wait_time)
                        lakeshore_temps.append(ctrl.read_temperature())
                        with open(filename1, 'a') as f1:
                            f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                                 str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                 str(format(float(lakeshore_temps[-1]), '.5f')) + '\t' +
                                 str(format(float(sc.get_cap()), '.5f')) + '\n')
                        if len(lakeshore_temps) > 120:
                            if (np.std(np.asarray(lakeshore_temps[-120:])) < lakeshore_stability) or (len(lakeshore_temps) > 7200):
                                print(f'Stabilized Lakeshore temperature at {ctrl.read_temperature()} K')
                                print(f'Lakeshore noise: {np.std(np.asarray(lakeshore_temps[-120:]))}')
                                caps = []
                                while True:
                                    time.sleep(wait_time)
                                    caps.append(sc.get_cap())
                                    with open(filename1, 'a') as f1:
                                        f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                                             str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                             str(format(float(ctrl.read_temperature()), '.5f')) + '\t' +
                                             str(format(float(caps[-1]), '.5f')) + '\n')
                                    if len(caps)>300: # mininum soak time of 5 minutes
                                        if (np.std(np.asarray(caps[-300:])) < cap_stability) or (len(caps) > 7200):
                                            print('STD: '+str(np.std(np.asarray(caps[-300:]))))
                                            if (np.std(np.asarray(caps[-300:])) < cap_stability):
                                                print('Stdev below accepted value')
                                            elif len(caps) > 7200: # maximum soak time of 2 hours
                                                print('Exceeded maximum soak time')
                                            print(f'Stabilized capacitance measurement, writing to file')
                                            with open(filename2, 'a') as f2:
                                                f2.write(str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                                     str(format(float(ctrl.read_temperature()), '.5f')) + '\t' +
                                                     str(format(np.mean(caps[-15:]), '.5f')) + '\n')
                                            break
                                        elif (np.std(np.asarray(caps[-300:])) < cap_stability*1.1) and (len(caps) > 3600):
                                            print('STD: '+str(np.std(np.asarray(caps[-300:]))))
                                            print('Stdev within 10 percent of accepted value after 1 hour')
                                            print(f'Stabilized capacitance measurement, writing to file')
                                            with open(filename2, 'a') as f2:
                                                f2.write(str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                                     str(format(float(ctrl.read_temperature()), '.5f')) + '\t' +
                                                     str(format(np.mean(caps[-15:]), '.5f')) + '\n')
                                            break
                                break
                    break

    cryo.warmup()
    print('Warming up cryostat')

def zero_strain_cell(sc, slew_rate=1, target_voltage=120, tol=0.1):
    '''
    carries out strain cell zeroing procedure as laid out in Razorbill documentatin. Energise both inner and outer (channel 1 and 2) stacks to +120V at room temperature with NO SAMPLE MOUNTED, and then allow stacks to return slowly to 0V by setting votlage to 0V and slew rate to 0.1, then turning off the outputs.

    args:
        - sc:           StrainClient object

    returns:
        - cap:          0 strain measurement of capacitance.
    '''
    sc.set_slew_rate(slew_rate)
    sc.set_output(1,1)
    sc.set_output(2,1)
    sc.set_voltage(1, target_voltage)
    sc.set_voltage(2, target_voltage)
    cond = True
    while cond:
        time.sleep(0.1)
        v1 = sc.get_voltage(1)
        v2 = sc.get_voltage(2)
        if (v1 > target_voltage-tol and v1 < target_voltage+tol) and (v2 > target_voltage-tol and v2 < target_voltage+tol):
            print('Reached 120 V on both channels')
            cond = False
    time.sleep(4)
    sc.set_slew_rate(0.1)
    sc.set_voltage(1,0)
    sc.set_voltage(2,0)
    cond = True
    while cond:
        time.sleep(0.1)
        v1 = sc.get_voltage(1)
        v2 = sc.get_voltage(2)
        if (v1 > -tol and v1 < tol) and (v2 > -tol and v2 < tol):
            cond=False
    sc.set_output(1,0)
    sc.set_output(2,0)
    cap = 0
    for i in range(100):
        cap = cap + sc.get_cap()
        time.sleep(0.1)
    return cap/100

######################
### Helper Methods ###
######################

def get_motor_range(start, end, step_size):
    '''
    helper function for returning a numpy array of positions corresponding to input. The convention is that step_size is always a positive number and that the motor moves from start to end.
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

    motors = list(map_dict.keys())

    # check validity and setup motor ranges and kwargs
    mranges = []
    mkwargs_dict = {}
    valid_motors = list(motor_dict.keys())
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

def initialize_motors(motors):
    mobj_dict = {}
    for m in motors:
        init_func = motor_dict[m]['init']
        mobj_dict[m] = init_func()
    return mobj_dict

def initialize_motors_threaded(motors):
    mobj_dict_locked = LockedDict({})
    initialization_threads = []
    for m in motors:
        init_func = motor_dict[m]['init']
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
    motors = list(mobj_dict.keys())
    for m in motors:
        obj = mobj_dict[m]
        close_func = motor_dict[m]['close']
        close_func(obj)

def move_motors_to_start(mobj_dict, mkwargs_dict, positions, print_flag=False):
    motors = list(mobj_dict.keys())
    for ii, m in enumerate(motors):
        move_back = motor_dict[m]['move_back']
        p = positions[0][ii]
        move_func = motor_dict[m]['move']
        read_func = motor_dict[m]['read']
        obj = mobj_dict[m]
        kwargs = mkwargs_dict[m]
        move_func(p - move_back, obj, **kwargs)
        move_func(p, obj, **kwargs)
        p_true = read_func(obj)
        if print_flag==True:
            print(f'Moved motor {m} to {p_true}.')

def move_motors(mobj_dict, mkwargs_dict, current_pos, start_pos, new_pos, print_flag=False):
    motors = list(mobj_dict.keys())
    for ii, m in enumerate(motors):
        p_old = current_pos[ii]
        p_new = new_pos[ii]
        p_start = start_pos[ii]
        if p_new!=p_old:
            move_func = motor_dict[m]['move']
            read_func = motor_dict[m]['read']
            obj = mobj_dict[m]
            kwargs = mkwargs_dict[m]
            if p_new==p_start:
                move_back = motor_dict[m]['move_back']
                move_func(p_new-move_back, obj, **kwargs)
                move_func(p_new, obj, **kwargs)
            else:
                move_func(p_new, obj, **kwargs)
            if print_flag==True:
                p_true = read_func(obj)
                print(f'Moved motor {m} to {p_true}.')

def read_motors(mobj_dict):
    motors = list(mobj_dict.keys())
    pos_dict = {}
    for m in motors:
        mobj = mobj_dict[m]
        read_func = motor_dict[m]['read']
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

def generate_filename(filename_head, filename, current_func):
    '''
    Generates filename_head and filename. If we supply it with no filename head or filename it autogenerates both. If we only supply a filaname it inferes the head to be current working directory.


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

def generate_metadata(mobj_dict):
    '''
    helper function that returns metadata dictionary
    '''

    metadata = {}
    for m in meta_motors:
        name = motor_dict[m]['name']
        try: # if in mobj_dict, don't close
            obj = mobj_dict[m]
            read = motor_dict[m]['read']
            pos = read(obj)
            metadata[name] = pos
        except: # if not in mobj_dict for whatever reason, try to initialize and then close
            try:
                obj = motor_dict[m]['init']()
                close = motor_dict[m]['close']
                read = motor_dict[m]['read']
                pos = read(obj)
                close(obj)
                metadata[name] = pos
            except:
                metadata[name]=None
    return metadata

def write_file_header(fname, header, metadata=None):
    '''
    Helper function that takes in metadata (a dictionary) and header (a list of strings) and writes data in the form:

    [Metadata]
    key: value
    [Data]
    header

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

#######################
### General Methods ###
#######################

def mapping(map_dict, single_point_function, function_args_dict):
    '''
    General purpose mapping function that evaluates single_point_function(args) at each point in the map. The map_dict has entries of the form 'axis':(start, end, num_steps, kwargs), where kwargs is a dictionary of key/value pairs appropriate for each motor 'move' function. For example, a temperature map might take the following map dictionary:

    map_dict = {'temp':(10,20,10,{'tolerance':0.01, 'wait_time':30})}

    UNDER CONSTRUCTION

    '''
    # Lock-in Amplifier initialization
    daq_objs = instrument_dict['zurich_lockin']['init']()

    # initialize axes
    axis_1 = motor_dict['axis_1']['init']()
    axis_2 = motor_dict['axis_1']['init']()

    # capture motor information and check for validity
    motors = list(map_dict.keys())
    for m in motors:
        valid_motors = list(motor_dict.keys())
        if m not in valid_motors:
            raise ValueError(f'Invalid motor name. Please select motors from the list {valid_motors}.')

    # initialize motors
    mobj_dict = {}
    for m in motors:
        init_func = motor_dict[m]['init']
        mobj_dict[m] = init_func()

    # setup motor ranges and kwargs - DO THIS RIGHT!
    mranges = []
    mkwargs_dict = {}
    for m in motors:
        start = map_dict[m][0]
        end = map_dict[m][1]
        nstep = map_dict[m][2]
        kwargs = map_dict[m][3]
        range = np.linspace(start, end, nstep)
        mranges.append(range)
        mkwargs_dict[m] = kwargs

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)

    # move motors to start position, using move_back to handle initial case
    for ii, m in enumerate(motors):
        move_back = motor_dict[m]['move_back']
        p = positions[0][ii] - move_back
        move_func = motor_dict[m]['move']
        obj = mobj_dict[m]
        kwargs = mkwargs_dict[m]
        move_func(p, obj, **kwargs)
        print(f'Moved motor {m} to {p}.')

    # loop over positions, only moving a motor if its target position has changed.
    current_pos = positions[0]
    for pos in positions:

            # move motors if position has changed
            for ii, m in enumerate(motors):
                p_old = current_pos[ii]
                p_new = pos[ii]
                if p_new!=p_old:
                    move_func = motor_dict[m][0]
                    obj = mobj_dict[m]
                    kwargs = mkwargs_dict[m]
                    move_func(p_new, obj, **kwargs)
                    print(f'Moved motor {m} to {p_new}.')

            # setup each filename - DO THIS RIGHT!
            totfilename = f'{filename_head}\{filename}_x{x_pos}_y{y_pos}.dat'

            # scan
            corotate_scan(num_steps, start_angle, end_angle, angle_offset, filename_head=filename_head, filename=totfilename, time_constant=time_constant, showplot=False, channel_index=channel_index, R_channel_index=R_channel_index, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2)

            current_pos = pos

    # close motors
    for m in motors:
        obj = mobj_dict[m]
        close_func = motor_dict[m]['close']
        close_func(obj)


#############
### Other ###
#############

class LockedVar:
    '''
    Minimal class to implement a locking variable. Contains two private attributes, a value and a lock, and a few methods for safetly reading writing value via the lock.
    '''
    def __init__(self, val):
        self._value = val
        self._lock = threading.Lock()

    def locked_read(self):
        with self._lock:
            return self._value

    def locked_update(self, val):
        with self._lock:
            self._value = val
class LockedDict:
    '''
    Minimal class to implement a locking dictionary. Contains two private attributes, a value and a lock, and a few methods for safetly reading writing value via the lock.
    '''
    def __init__(self, dict):
        self._dict = dict
        self._lock = threading.Lock()

    def locked_read(self, key):
        with self._lock:
            return self._dict[key]

    def locked_update(self, key, val):
        with self._lock:
            self._dict[key] = val

    def locked_get_dict(self):
        with self._lock:
            return self._dict
