'''
Main file for controlling lab equipment and orchestrating measurements, with a specific eye to procedures
'''

from strain_control.strain_client import StrainClient
import OrensteinLab_git.Measurement.Alex.control as ctrl
import OrensteinLab_git.Instrument.montana.cryocore as cryocore
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import threading
from tqdm.auto import tqdm
'''
Features to add:
    - plotting for motor_scans
    - generalize handling of writing files (headers, metadata, etc) using motor names and other such features that can be easily set at the top of the file.
    - lockin_time_series can read motors and add to the file/plots
    - add strain cell as a motor - needs some thought
    - default kwargs for each motor when running from meas.
    - record only measured values, not setpoints

'''

'''
map_dict are dictionaries with entries of the form {'motor':(start, stop, step, kwargs)}, where kwargs is another dictionary with form {'key':value} where key is the name of kwarg and value the desired value.

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
lockin_header = ['Demod x', 'Demod y', 'r', 'Demod x_R', 'Demod y_R', 'Demod r_R']

#####################
### System Motors ###
#####################
# entries of the form motor:{'move':move_function, 'read':read_function, 'init':initialize_function, 'close':close_function}
#'coil':(ctrl.set_coil, ctrl.initialize_coil, ctrl.close_coil),
# consider changing to yaml file, can also be used for easily have default configurations for instruments.
motor_dict = {
'x':{'move':ctrl.move_x, 'read':ctrl.read_x, 'init':ctrl.initialize_attocube, 'close':ctrl.close_attocube, 'move_back':10, 'name':'x (um)'},

'y':{'move':ctrl.move_y, 'read':ctrl.read_y, 'init':ctrl.initialize_attocube, 'close':ctrl.close_attocube, 'move_back':10, 'name':'y (um)'},

'z':{'move':ctrl.move_z, 'read':ctrl.read_z, 'init':ctrl.initialize_attocube, 'close':ctrl.close_attocube, 'move_back':10, 'name':'z (um)'},

'temp':{'move':ctrl.set_temperature, 'read':ctrl.read_temperature, 'init':ctrl.initialize_lakeshore, 'close':ctrl.close_lakeshore, 'move_back':0, 'name':'Temperature (K)'},

'axis_1':{'move':ctrl.rotate_axis_1, 'read':ctrl.read_axis_1, 'init':ctrl.initialize_rot_axis_1, 'close':ctrl.close_rot_axis_1, 'move_back':1, 'name':'Angle 1 (deg)'},

'axis_2':{'move':ctrl.rotate_axis_2, 'read':ctrl.read_axis_2, 'init':ctrl.initialize_rot_axis_2, 'close':ctrl.close_rot_axis_2, 'move_back':1, 'name':'Angle 2 (deg)'},

'delay_stage':{'move':ctrl.move_delay_stage, 'read':ctrl.read_delay_stage, 'init':ctrl.initialize_delay_stage, 'close':ctrl.close_delay_stage, 'move_back':1, 'name':'Delay Stage (mm)'},

'strain_cap':{'move':ctrl.set_strain_capacitance, 'read':ctrl.read_strain_capacitance, 'init':ctrl.initialize_strain_cell_client, 'close':ctrl.close_strain_cell_client, 'move_back':0, 'name':'Capacitance (pF)'},

'strain_ps':{'move':ctrl.set_strain_ps, 'read':ctrl.read_strain_ps, 'init':ctrl.initialize_strain_cell_client, 'close':ctrl.close_strain_cell_client, 'move_back':0, 'name':'Voltage (V)'},

'zurich_output':{'move':ctrl.set_zurich_output_amplitude, 'read':ctrl.read_zurich_output_amplitude, 'init':ctrl.initialize_zurich_lockin, 'close':ctrl.close_zurich_lockin, 'move_back':0, 'name':'Lock-in Ouput Voltage (V)'},

'zurich_frequency':{'move':ctrl.set_zurich_frequency, 'read':ctrl.read_zurich_frequency, 'init':ctrl.initialize_zurich_lockin, 'close':ctrl.close_zurich_lockin, 'move_back':0, 'name':'Lock-in Frequency (Hz)'}
}

##########################
### System Instruments ###
##########################
instrument_dict = {
'zurich_lockin':{'read':ctrl.read_zurich_lockin, 'init':ctrl.initialize_zurich_lockin}
}

################
### Medthods ###
################

def lockin_time_series(recording_time, filename_head=None, filename=None, measure_motors=[], time_constant=0.3, channel_index=1, R_channel_index=2):
    '''
    aquires data on the lockin over a specified length of time.
    '''
    # initialize zurich lockin and setup read function
    daq_objs  = instrument_dict['zurich_lockin']['init']()
    read_lockin = instrument_dict['zurich_lockin']['read']

    # initialize additional motors to measure during scan
    mobj_dict = initialize_motors(measure_motors)

    # initialize data bins
    time_record = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # setup plot
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
    if filename_head!=None and filename!=None:
        fname = get_unique_filename(filename_head, filename)
        header_motors = [motor_dict[m]['name'] for m in measure_motors]
        header = ['Time (s)', 'Demod x', 'Demod y', 'R']+header_motors
        write_file_header(fname, header)

    # loop
    t_delay = 0
    tic = time.perf_counter()
    while (t_delay<recording_time):
        time.sleep(time_constant*4)
        toc = time.perf_counter()
        t_delay = toc - tic

        x, y, r, x_R, y_R, r_R = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index)

        # read additional motors
        measured_positions_dict = read_motors(mobj_dict)
        measured_positions = [measured_positions_dict[m] for m in measure_motors]

        time_record = np.append(time_record, t_delay)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # update plot
        draw_x.set_data(time_record-time_record[0],demod_x)
        draw_y.set_data(time_record-time_record[0],demod_y)
        draw_r.set_data(time_record-time_record[0],demod_r)
        for ax in axes:
            ax.relim()
            ax.autoscale()
        fig.canvas.draw()
        fig.canvas.flush_events()

        # update file
        if filename_head!=None and filename!=None:
            vars = [t_delay, x, y, r]+measured_positions
            append_data_to_file(fname, vars)

    return time_record, demod_x, demod_y, demod_r

def rotate_scan(start_angle, end_angle, step_size, filename_head=None, filename=None, axis_index=1, measure_motors=[], showplot=True, time_constant=0.3, channel_index=1, R_channel_index=2, daq_objs=None, axis_1=None, axis_2=None):

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
    mobj_dict = initialize_motors(measure_motors)

    # setup measureables
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # convert input to angle lists
    angles = get_motor_range(start_angle, end_angle, step_size)

    # setup measureables
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # initialize file
    if filename_head!=None and filename!=None:
        fname = get_unique_filename(filename_head, filename)
        header_motors = [motor_dict[m]['name'] for m in measure_motors]
        header = [motor_dict['axis_1']['name'], motor_dict['axis_2']['name']]+lockin_header+header_motors
        write_file_header(fname, header)

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
        x, y, r, x_R, y_R, r_R = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index)
        angle_pos_1 = read_axis_1(axis=axis_1, print_flag=False)
        angle_pos_2 = read_axis_2(axis=axis_2, print_flag=False)
        if axis_index==1:
            angle_pos = angle_pos_1
        elif axis_index==2:
            angle_pos = angle_pos_2

        # read additional motors
        measured_positions_dict = read_motors(mobj_dict)
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
        if filename_head!=None and filename!=None:
            vars = [angle_pos_1, angle_pos_2, x, y, r, x_R, y_R, r_R]+measured_positions
            append_data_to_file(fname, vars)

    # move motors back to original positions
    move_axis(start_angle, axis=axis)
    #move_other_axis(start_angle, axis=other_axis)

    return position, demod_x, demod_y, demod_r

def corotate_scan(start_angle, end_angle, step_size, angle_offset, rate_axis_2=1, filename_head=None, filename=None, measure_motors=[], showplot=True, time_constant=0.3, channel_index=1, R_channel_index=2, daq_objs=None, axis_1=None, axis_2=None):
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
    mobj_dict = initialize_motors(measure_motors)

    # setup measureables
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # convert input to angle lists
    angles_1 = get_motor_range(start_angle, end_angle, step_size)
    angles_2 = rate_axis_2*angles_1 + angle_offset

    # setup measureables
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    # initialize file
    if filename_head!=None and filename!=None:
        fname = get_unique_filename(filename_head, filename)
        header_motors = [motor_dict[m]['name'] for m in measure_motors]
        header = header = [motor_dict['axis_1']['name'], motor_dict['axis_2']['name']]+lockin_header+header_motors
        write_file_header(fname, header)

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
            ctrl.corotate_axes(1, 2, angle_1-move_back_1, angle_2-move_back_2, axis_1=axis_1, axis_2=axis_2)
            #ctrl.corotate_axes(1, 2, angle_1, angle_2, axis_1=axis_1, axis_2=axis_2)
            time.sleep(2) # was set to 2 but I don't think this matters much
        ctrl.corotate_axes(1, 2, angle_1, angle_2, axis_1=axis_1, axis_2=axis_2)
        #move_axis_1(angle_1, axis=axis_1)
        #move_axis_2(angle_2, axis=axis_2)
        time.sleep(time_constant)

        # read lockin, rotators
        x, y, r, x_R, y_R, r_R = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index)
        angle_pos_1 = read_axis_1(axis=axis_1, print_flag=False)
        angle_pos_2 = read_axis_2(axis=axis_2, print_flag=False)

        # read additional motors
        measured_positions_dict = read_motors(mobj_dict)
        measured_positions = [measured_positions_dict[m] for m in measure_motors]

        position = np.append(position, angle_pos_1)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # write to file
        if filename_head!=None and filename!=None:
            vars = [angle_pos_1, angle_pos_2, x, y, r, x_R, y_R, r_R]+measured_positions
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
    ctrl.corotate_axes(1, 2, start_angle, start_angle+angle_offset, axis_1=axis_1, axis_2=axis_2)

    return position, demod_x, demod_y, demod_r

def motor_scan(map_dict, filename_head=None, filename=None, measure_motors=[], showplot=True, time_constant=0.3, channel_index=1, R_channel_index=3, print_flag=False):
    '''
    utility to record lockin measurement as a function of motors specified by dictionary map_dict.
    '''

    # Lock-in Amplifier initialization
    daq_objs = instrument_dict['zurich_lockin']['init']()
    read_lockin = instrument_dict['zurich_lockin']['read']

    # capture motor information and initialize
    motors, mranges, mkwargs_dict = capture_motor_information(map_dict)
    mobj_dict = initialize_motors(motors)
    mobj_measure_dict = initialize_motors(measure_motors)

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)
    #print(positions)

    # setup file with header
    if filename_head!=None and filename!=None:
        fname = get_unique_filename(filename_head, filename)
        header_motors = [motor_dict[m]['name'] for m in measure_motors]
        header = [motor_dict[m]['name'] for m in motors]+[motor_dict[m]['name']+str(' measured') for m in motors]+lockin_header+header_motors
        write_file_header(fname, header)

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

        # acquire data
        lockin_meas = read_lockin(daq_objs=daq_objs, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index)
        x, y, r, x_R, y_R, r_R = lockin_meas

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
        if filename_head!=None and filename!=None:
            append_data_to_file(fname, list(pos)+list(measured_positions)+list(lockin_meas)+list(measured_motors_positions))

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

        current_pos = pos

    # close motors
    close_motors(mobj_dict)

def rotate_map(map_dict, start_angle, end_angle, step_size, filename_head=None, filename=None, axis_index=1, measure_motors=[], showplot=False, time_constant=0.3, channel_index=1, R_channel_index=2, daq_objs=None, print_flag=False):

    # Lock-in Amplifier initialization
    daq_objs = instrument_dict['zurich_lockin']['init']()

    # initialize rotation axes
    axis_1 = motor_dict['axis_1']['init']()
    axis_2 = motor_dict['axis_2']['init']()

    # capture motor information and initialize
    motors, mranges, mkwargs_dict = capture_motor_information(map_dict)
    mobj_dict = initialize_motors(motors)

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)

    # move motors to start position, using move_back to handle initial case
    move_motors_to_start(mobj_dict, mkwargs_dict, positions, print_flag=print_flag)

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
        rotate_scan(start_angle, end_angle, step_size, filename_head=filename_head, filename=expanded_filename, axis_index=axis_index, measure_motors=measure_motors, showplot=showplot, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2)

        current_pos = pos

    # close motors
    close_motors(mobj_dict)

def corotate_map(map_dict, start_angle, end_angle, step_size, angle_offset, rate_axis_2=1, filename_head=None, filename=None, measure_motors=[], showplot=False, time_constant=0.3, channel_index=1, R_channel_index=2, print_flag=False):
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
    mobj_dict = initialize_motors(motors)

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)

    # move motors to start position, using move_back to handle initial case
    move_motors_to_start(mobj_dict, mkwargs_dict, positions, print_flag=print_flag)

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
        corotate_scan(start_angle, end_angle, step_size, angle_offset, rate_axis_2=rate_axis_2, filename_head=filename_head, filename=expanded_filename, measure_motors=measure_motors, showplot=showplot, time_constant=time_constant, channel_index=channel_index, R_channel_index=R_channel_index, daq_objs=daq_objs, axis_1=axis_1, axis_2=axis_2)

        current_pos = pos

    # close motors
    close_motors(mobj_dict)

def find_balance_angle(start_angle, end_angle, step_size, balance_at=0, go_to_balance_angle=True, axis_index=2):
    '''
    Assuming we are measuring in DC mode above a transition or on GaAs, carries out a rotate_scan. Find angle by carrying out a linear fit, such that the angle range should be taken to be very small.

    By default, moves axis 2 to find balance angle.

    If go_to_balance_angle is set to true, moves stages to
    '''

    # initialize axes
    axis_1 = motor_dict['axis_1']['init']()
    axis_2 = motor_dict['axis_2']['init']()
    move_axis_1 = motor_dict['axis_1']['move']
    move_axis_2 = motor_dict['axis_2']['move']
    move_back_1 = motor_dict['axis_1']['move_back']
    move_back_2 = motor_dict['axis_2']['move_back']

    # move both motors to 0
    move_axis_1(balance_at-move_back_1)
    move_axis_2(balance_at-move_back_2)
    move_axis_1(balance_at)
    move_axis_2(balance_at)

    positions, demod_x, demod_y, demod_r = rotate_scan(balance_at+start_angle, balance_at+end_angle, step_size, axis_index=axis_index)

    # linear fit
    fit_params = np.polyfit(positions, demod_x, 1)
    balance_angle = -fit_params[1]/fit_params[0]
    angles_vect = np.linspace(start_angle, end_angle, 1000)
    fit = fit_params[0]*angles_vect + fit_params[1]

    # display result
    fig, ax = plt.subplots(1)
    ax.set_ylabel('Angle (deg)')
    ax.set_ylabel('Demod X')
    ax.plot(positions, demod_x, 'o', ms=5, color='blue')
    ax.plot(angles_vect, fit, '-', color='black')

    # set balance angle relative to 0
    balance_angle = balance_angle - balance_at

    if go_to_balance_angle: # is there a good way to automatically re-zero zero?
        if axis_index==1:
            move_axis_1(balance_angle)
        else:
            move_axis_2(balance_angle)

    print(f'Balance angle: {balance_angle}')
    return balance_angle

def align_delay_stage(wait_time=5):

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
            move(-125, delay_stage)
            time.sleep(wait)
            move(125, delay_stage)
            with lock: # read locked variable
                run_cond = run
                # print(run_cond)

    # start thread
    move_ds_thread = threading.Thread(target=osc_delay_stage, args=(wait_time,))
    move_ds_thread.start()

    # ask to stop
    answer = input('Stop?')

    # shutdown
    with lock:
        run = False
    move_ds_thread.join()

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
            p_true = read_func(obj)
            if print_flag==True:
                print(f'Moved motor {m} to {p_true}.')

def read_motors(mobj_dict):
    motors = list(mobj_dict.keys())
    pos_dict = {}
    for ii, m in enumerate(motors):
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

def get_unique_filename(filename_head, filename):
    fname = f'{filename_head}\\{filename}.dat'
    if not os.path.exists(fname):
        return fname

    path, name = os.path.split(fname)
    name, ext = os.path.splitext(name)

    make_fn = lambda i: os.path.join(path, '%s_%s%s' % (name, i, ext))

    # dir
    # for fname in dir:
    #     if filename==fname:
    #         match = re.match(filename, fname)
    #         int(bool(match))

    for i in range(1, 1000):
        uni_fn = make_fn(i)
        if not os.path.exists(uni_fn):
            return uni_fn

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
            for key in metadata:
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

def get_all_motor_positions():
    '''
    helper function that reads all positions of motors in motor_dict. For use in making metada.
    '''
    return 1

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
