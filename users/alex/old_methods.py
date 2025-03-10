'''
This file contains old balancing methods that I wrote at some point but which I never really used.
'''

def motor_scan_balance(map_dict, balance, balance_table=None, slope=0, tol=0, balance_channel=3, autobalance_flag=True, filename_head=None, filename=None, showplot=True, time_constant=0.3, channel_index=1, R_channel_index=4, print_flag=False, savefile=True, metadata={}, daq_objs=None):
    '''
    utility to record lockin measurement as a function of motors specified by dictionary map_dict.

    autobalance before each scan of last motor in map_dict according to balance table, which has the form of an xarray object with balance angle data over the various coordinates.

    only works if corotate_axes12 is in the map_dict (not sure that this is true anymore...)
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
        elif 'axis_1' in list(mobj_measure_dict.keys()):
            axis_1 = mobj_measure_dict['axis_1']
        else:
            axis_1 = motor_dict['axis_1']['init']()
        if 'axis_2' in motors:
            axis_2 = mobj_dict['axis_2']
        elif 'axis_2' in list(mobj_measure_dict.keys()):
            axis_2 = mobj_measure_dict['axis_2']
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

def find_balance_angle_macro(start_angle, end_angle, step_size, step_size_fine=0.2, window=3, axis_index=2, balance_at=0, channel_index=1, time_constant=0.3, daq_objs=None, axis_1=None, axis_2=None):

    # initialize zurich lockin and setup read function
    if daq_objs==None:
        init_func = INSTRUMENT_DICT['zurich_lockin']['init']
        daq_objs = init_func()
    read_lockin = INSTRUMENT_DICT['zurich_lockin']['read']

    # initialize axes and setup move functions
    if axis_1 == None:
        init_func = MOTOR_DICT['axis_1']['init']
        axis_1 = init_func()
    if axis_2 == None:
        init_func = MOTOR_DICT['axis_2']['init']
        axis_2 = init_func()

    # initialize axes funcs
    move_axis_1 = MOTOR_DICT['axis_1']['move']
    move_axis_2 = MOTOR_DICT['axis_2']['move']
    move_back_1 = MOTOR_DICT['axis_1']['move_back']
    move_back_2 = MOTOR_DICT['axis_2']['move_back']

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

def autobalance_cont(slope, tolerance=None, daq_objs=None, axis_1=None, axis_2=None, balance_at=None, offset=0, channel_index=1, time_constant=0.3, print_flag=True, lockin_index=0):

    # initialize lockin
    if daq_objs==None:
        daq_objs = INSTRUMENT_DICT['zurich_lockin']['init']()
    else:
        pass
    read_lockin = INSTRUMENT_DICT['zurich_lockin']['read']

    # initialize axes and setup axis move functions
    if axis_1 == None:
        init_func = MOTOR_DICT['axis_1']['init']
        axis_1 = init_func()
    if axis_2 == None:
        init_func = MOTOR_DICT['axis_2']['init']
        axis_2 = init_func()
    move_axis_1 = MOTOR_DICT['axis_1']['move']
    move_axis_2 = MOTOR_DICT['axis_2']['move']
    read_axis_1 = MOTOR_DICT['axis_1']['read']
    read_axis_2 = MOTOR_DICT['axis_2']['read']
    move_back_1 = MOTOR_DICT['axis_1']['move_back']
    move_back_2 = MOTOR_DICT['axis_2']['move_back']

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
        axis_1 = MOTOR_DICT['axis_1']['init']()
    if 'axis_2' in motors:
        axis_2 = mobj_dict['axis_2']
    else:
        axis_2 = MOTOR_DICT['axis_2']['init']()

    # generate positions recursively
    positions = gen_positions_recurse(mranges, len(mranges)-1)
    #print(positions)

    # setup file for writing
    filename_head, filename = generate_filename(filename_head, filename, 'balance_table')
    for m in motors:
        filename = filename+f'_{m}'
    fname = get_unique_filename(filename_head, filename)
    header = [MOTOR_DICT[m]['name'] for m in motors]+['Balance Angle (deg)', 'Slope', 'Tolerance']
    write_file_header(fname, header, metadata)

    # move motors to start position, using move_back to handle initial case
    move_motors_to_start(mobj_dict, mkwargs_dict, positions, print_flag=print_flag)

    # loop over positions, only moving a motor if its target position has changed.
    if 'axis_1' in motors:
        angle1_idx = motors.index('axis_1')
    else:
        balance_at = MOTOR_DICT['axis_1']['read']()
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
    motor_names = [MOTOR_DICT[m]['name'] for m in motors]
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

    daq_objs = INSTRUMENT_DICT['zurich_lockin']['init']()

    # initialize rotation axes
    axis_1 = MOTOR_DICT['axis_1']['init']()
    axis_2 = MOTOR_DICT['axis_2']['init']()

    # setup file for writing
    filename_head, filename = generate_filename(filename_head, filename, 'balance_table')
    fname = get_unique_filename(filename_head, filename)
    header = ['Corotation Axes (deg)']+['Balance Angle (deg)', 'Slope', 'Tolerance']
    write_file_header(fname, header, metadata)

    positions = get_motor_range(start_angle, end_angle, step_size)

    # move motors to start position
    MOTOR_DICT['axis_1']['move'](0, axis_1)
    MOTOR_DICT['axis_2']['move'](0+bal_angle_init, axis_2)

    bal_angle = bal_angle_init
    for ii in tqdm(range(len(positions))):
        pos = positions[ii]

        # move motors
        MOTOR_DICT['axis_1']['move'](pos, axis_1)
        MOTOR_DICT['axis_2']['move'](pos+bal_angle, axis_2)

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


###
### Helper methods
###

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