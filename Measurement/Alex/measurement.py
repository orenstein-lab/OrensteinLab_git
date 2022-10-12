'''
Main file for controlling lab equipment and orchestrating measurements, with a specific eye to procedures
'''

from strain_control.strain_client import StrainClient
import OrensteinLab_git.Measurement.Alex.control as ctrl
import OrensteinLab_git.Instrument.montana.cryocore as cryocore
import time
import numpy as np
import matplotlib.pyplot as plt

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
            for key, value in metadata.items():
                f.write(f'{key}:\t{value}\n')
            f.write('[DATA]\n')
        for item in header:
            f.write(str(item)+'\t')
        f.write('\n')
        for line in data:
            for item in line:
                f.write(str(item)+'\t')
            f.write('\n')

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

'''
Features to add:
    - quick motor scanning functionality - ideally one tool
    - measure power by averaging over lock-in.
    - eliminate _quick functions - just make them one function with some flexible flags.
    - change how points are chosen
    - make code more readible.
'''

def Corotate_quick(num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, step_size_2, go_back_2, channel_index, time_constant):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    end_pos_1 = start_pos_1 + step_size_1 * (num_of_steps-1)
    end_pos_2 = start_pos_2 + step_size_2 * (num_of_steps-1)
    scan_range_1 = np.arange(start_pos_1, start_pos_1 + step_size_1 * num_of_steps, step_size_1)
    scan_range_2 = np.arange(start_pos_2, start_pos_2 + step_size_2 * num_of_steps, step_size_2)
    axis_rot_1 = newport.NewportESP301Axis(controller,axis_index_1-1)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller,axis_index_2-1)
    axis_rot_2.enable()

    global button
    button = True

    #Quantities to be measured
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    fig = plt.figure(figsize=(8,10))
    gs = fig.add_gridspec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax1.set_xlabel('Angle (deg)')
    ax1.set_ylabel('Demod x')
    ax2.set_xlabel('Angle (deg)')
    ax2.set_ylabel('Demod y')
    ax3.set_xlabel('Angle (deg)')
    ax3.set_ylabel('R')
    draw_x, = ax1.plot([],'-o')
    draw_y, = ax2.plot([],'-o')
    draw_r, = ax3.plot([],'-o')
    fig.canvas.draw()
    fig.show()

    while True:
        time.sleep(0.03)
        try:
            if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                break
        except ValueError:
            pass
    #Scan
    for i in np.arange(0,len(scan_range_1),1):
        if button == False:
            break
        pos_1 = scan_range_1[i]
        pos_2 = scan_range_2[i]
        if (pos_1 == start_pos_1):
            axis_rot_1.move(pos_1-go_back_1,absolute=True)
            axis_rot_2.move(pos_2-go_back_2,absolute=True)
            while True:
                time.sleep(0.03)
                try:
                    if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                        break
                except ValueError:
                    pass
        axis_rot_1.move(pos_1,absolute=True)
        time.sleep(0.03)
        axis_rot_2.move(pos_2,absolute=True)
        while True:
            time.sleep(0.03)
            try:
                if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                    break
            except ValueError:
                pass
        time.sleep(time_constant*4)
        sample = daq.getSample(channel_name[channel_index-1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        position = np.append(position, axis_rot_1.position)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        #Plot
        draw_x.set_data(position,demod_x)
        draw_y.set_data(position,demod_y)
        draw_r.set_data(position,demod_r)
        ax1.relim()
        ax1.autoscale()
        ax2.relim()
        ax2.autoscale()
        ax3.relim()
        ax3.autoscale()
        fig.canvas.draw()
        fig.canvas.flush_events()

    axis_rot_1.move(start_pos_1,absolute=True)
    axis_rot_2.move(start_pos_2,absolute=True)

def Corotate(num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, step_size_2, go_back_2, channel_index, R_channel_index, time_constant, showplot, filename_head, filename):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    totfilename = filename_head + '\\' + filename + '.dat'
    file = open(totfilename,'a')
    file.write("Angle_1 (deg)"+'\t'+"Angle_2 (deg)"+'\t'+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
    file.close()

    end_pos_1 = start_pos_1 + step_size_1 * (num_of_steps-1)
    end_pos_2 = start_pos_2 + step_size_2 * (num_of_steps-1)
    scan_range_1 = np.arange(start_pos_1, start_pos_1 + step_size_1 * num_of_steps, step_size_1)
    scan_range_2 = np.arange(start_pos_2, start_pos_2 + step_size_2 * num_of_steps, step_size_2)
    axis_rot_1 = newport.NewportESP301Axis(controller,axis_index_1-1)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller,axis_index_2-1)
    axis_rot_2.enable()

    global button
    button = True
    # thread1 = threading.Thread(target=get_input)
    # thread1.start()

    #Quantities to be measured
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    if showplot == True:
        fig = plt.figure(figsize=(8,10))
        gs = fig.add_gridspec(3, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[2, 0])
        ax1.grid(True)
        ax2.grid(True)
        ax3.grid(True)
        ax1.set_xlabel('Angle (deg)')
        ax1.set_ylabel('Demod x')
        ax2.set_xlabel('Angle (deg)')
        ax2.set_ylabel('Demod y')
        ax3.set_x

        label('Angle (deg)')
        ax3.set_ylabel('R')
        draw_x, = ax1.plot([],'-o')
        draw_y, = ax2.plot([],'-o')
        draw_r, = ax3.plot([],'-o')
        fig.canvas.draw()
        fig.show()

    #Scan
    while (axis_rot_1.is_motion_done == False) or (axis_rot_2.is_motion_done == False):
        time.sleep(0.1)

    for i in np.arange(0,len(scan_range_1),1):
        if button == False:
            break
        pos_1 = scan_range_1[i]
        pos_2 = scan_range_2[i]
        if (pos_1 == start_pos_1):
            axis_rot_1.move(pos_1-go_back_1,absolute=True)
            axis_rot_2.move(pos_2-go_back_2,absolute=True)
            while True:
                try:
                    axis_rot_1.is_motion_done or axis_rot_2.is_motion_done
                    break
                except ValueError:
                    pass
        axis_rot_1.move(pos_1,absolute=True)
        time.sleep(0.03)
        axis_rot_2.move(pos_2,absolute=True)
        while True:
            try:
                axis_rot_1.is_motion_done or axis_rot_2.is_motion_done
                break
            except ValueError:
                pass
        time.sleep(time_constant*4)
        sample = daq.getSample(channel_name[channel_index-1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        sample_R = daq.getSample(channel_name[R_channel_index - 1] % device)
        sample_R["R"] = np.abs(sample_R["x"] + 1j * sample_R["y"])
        x_R = sample_R["x"][0]
        y_R = sample_R["y"][0]
        r_R = sample_R["R"][0]

        position = np.append(position, axis_rot_1.position)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        #Plot
        if showplot == True:
            draw_x.set_data(position,demod_x)
            draw_y.set_data(position,demod_y)
            draw_r.set_data(position,demod_r)
            ax1.relim()
            ax1.autoscale()
            ax2.relim()
            ax2.autoscale()
            ax3.relim()
            ax3.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

        #Save file
        file = open(totfilename,'a')
        file.write(format(axis_rot_1.position,'.15f')+"\t"+format(axis_rot_2.position,'.15f')+"\t"+format(x,'.15f')+'\t'+format(y,'.15f')+'\t'+format(r,'.15f')+'\t'+format(x_R,'.15f')+'\t'+format(y_R,'.15f')+'\t'+format(r_R,'.15f')+'\n')
        file.close()
    #thread1.join()
    axis_rot_1.move(start_pos_1,absolute=True)
    axis_rot_2.move(start_pos_2,absolute=True)

def Corotate_map_scan(x_start, x_step, x_num, y_start, y_step, y_num, num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, go_back_2, channel_index, time_constant, filename_head, filename):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    #Attocube initialization
    ax = {'x':0,'y':1,'z':2}
    anc = Positioner()

    end_pos_1 = start_pos_1 + step_size_1 * (num_of_steps-1)
    scan_range_1 = np.arange(start_pos_1, start_pos_1 + step_size_1 * num_of_steps, step_size_1)
    #scan_range_2 = np.arange(start_pos_2, start_pos_2 + step_size_2 * num_of_steps, step_size_2)
    axis_rot_1 = newport.NewportESP301Axis(controller,axis_index_1-1)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller,axis_index_2-1)
    axis_rot_2.enable()
    start_pos = start_pos_1
    end_pos = start_pos_1 + step_size_1 * (num_of_steps-1)

    x_end = x_start + x_step * x_num
    y_end = y_start + y_step * y_num

    x_range = np.arange(x_start, x_end, x_step)
    y_range = np.arange(y_start, y_end, y_step)
    go_back = 10
    x_tor = 1
    y_tor = 1

    start_pos_2 = float(start_pos_2)
    scan_range_2 = np.arange(start_pos_2, start_pos_2 + step_size_1 * num_of_steps, step_size_1)
    end_pos_2 = start_pos_2 + step_size_1 * (num_of_steps-1)

    for y_pos in y_range:
        for x_pos in x_range:
            if (x_pos == x_range[0] and y_pos == y_range[0]):
                x_target = x_pos-go_back
                y_target = y_pos-go_back
                anc.moveAbsolute(ax['x'], int(x_target*1000))
                anc.moveAbsolute(ax['y'], int(y_target*1000))
                x_error = np.abs(x_target-anc.getPosition(ax['x'])/1000)
                y_error = np.abs(y_target-anc.getPosition(ax['y'])/1000)
                while (x_error >= x_tor) or (y_error >= y_tor):
                    time.sleep(0.1)
                    x_error = np.abs(x_target-anc.getPosition(ax['x'])/1000)
                    y_error = np.abs(y_target-anc.getPosition(ax['y'])/1000)
                    if (x_error >= x_tor):
                        anc.moveAbsolute(ax['x'], int(x_target*1000))
                    if (y_error >= y_tor):
                        anc.moveAbsolute(ax['y'], int(y_target*1000))
            anc.moveAbsolute(ax['x'], int(x_pos*1000))
            x_error = np.abs(x_pos-anc.getPosition(ax['x'])/1000)
            while (x_error >= x_tor):
                time.sleep(0.1)
                #clear_output(wait=True)
                x_error = np.abs(x_pos-anc.getPosition(ax['x'])/1000)
                #print(x_error)
                if (x_error >= x_tor):
                    anc.moveAbsolute(ax['x'], int(x_pos*1000))
            anc.moveAbsolute(ax['y'], int(y_pos*1000))
            y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
            while (y_error >= y_tor):
                time.sleep(0.1)
                #clear_output(wait=True)
                y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
                #print(y_error)
                if (y_error >= y_tor):
                    anc.moveAbsolute(ax['y'], int(y_pos*1000))
            print("moved to "+str(anc.getPosition(ax['x'])/1000)+","+str(anc.getPosition(ax['y'])/1000))


            totfilename = filename_head + '\\'+ filename + '_x' + str(x_pos) + '_y' + str(y_pos) +'.dat'
            file = open(totfilename,'a')
            file.write("Angle 1 (deg)"+'\t'+"Angle 2 (deg)"+'\t'+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
            file.close()

            #Quantities to be measured
            position = np.array([])
            demod_x = np.array([])
            demod_y = np.array([])
            demod_r = np.array([])

            #Scan
            for i in np.arange(0,len(scan_range_1),1):
                pos_1 = scan_range_1[i]
                pos_2 = scan_range_2[i]
                if (pos_1 == start_pos_1):
                    axis_rot_1.move(pos_1-go_back_1,absolute=True)
                    axis_rot_2.move(pos_2-go_back_2,absolute=True)
                    while True:
                        time.sleep(0.03)
                        try:
                            if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                                break
                        except ValueError:
                            pass
                axis_rot_1.move(pos_1,absolute=True)
                time.sleep(0.03)
                axis_rot_2.move(pos_2,absolute=True)
                while True:
                    time.sleep(0.03)
                    try:
                        if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                            break
                    except ValueError:
                        pass
                time.sleep(time_constant*4)
                sample = daq.getSample(channel_name[channel_index-1] % device)
                sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
                x = sample["x"][0]
                y = sample["y"][0]
                r = sample["R"][0]
                position = np.append(position, axis_rot_1.position)
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_r = np.append(demod_r, r)

                #Save file
                file = open(totfilename,'a')
                file.write(format(float(position[-1]), '.15f')+"\t"+format(float(position[-1]+start_pos_2), '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
                file.close()

    anc.close()

def Corotate_map_Tdep(x_start_pos, x_step_size, x_num_steps, y_start_pos, y_step_size, y_num_steps, num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, go_back_2, start_temp, end_temp, temp_step_size, channel_index, time_constant, filename_head, filename):
    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    # Lakeshore initialization
    import OrensteinLab.Instrument.Lakeshore.Lakeshore335 as ls
    lsobj = ls.initialization_lakeshore335()
    ls.set_ramp(lsobj, 1, 0, 0)

    # Set time constant as specified
    daq.setDouble('/%s/demods/0/timeconstant' % device, time_constant)

    # List of temperatures to measure
    Trange = get_temp_values(start_temp, end_temp, temp_step_size)

    x_start = x_start_pos
    x_step = x_step_size
    x_num = x_num_steps
    y_start = y_start_pos
    y_step = y_step_size
    y_num = y_num_steps
    # num_of_steps = num_of_steps


    for temp in Trange:
        # Change temperature
        ls.set_setpoint(lsobj, 1, temp)
        time.sleep(0.1)
        currT = []
        for m in range(60):
            currT.append(ls.read_temperature(lsobj))
            if m >= 2 and abs(np.mean(currT[-3:]) - temp) < 0.05:
                time.sleep(80)
                break
            else:
                time.sleep(1)
        print(ls.read_temperature(lsobj))

        showplot = False
        temp_filename = filename + '_' + str(temp)
        Corotate_map_scan(x_start, x_step, x_num, y_start, y_step, y_num, num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, go_back_2, channel_index, time_constant, filename_head, temp_filename)

    print("Scans finished!")
    ls.close_lakeshore335(lsobj)
