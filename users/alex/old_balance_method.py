'''
This file contains old balancing methods that I wrote at some point but which I never really used.
'''

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
