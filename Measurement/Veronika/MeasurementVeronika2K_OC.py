#!/usr/bin/env python
# coding: utf-8

# In[ ]:



get_ipython().run_line_magic('matplotlib', 'notebook')
import zhinst.utils as ziutils
import instruments.newport as newport
import OrensteinLab.Instrument.OptiCool.OptiCool_Control as optc
import OrensteinLab.Instrument.Lakeshore.Lakeshore335 as ls
import numpy as np
import os
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os.path
import time
import datetime
from IPython.display import clear_output
import threading
import ipywidgets as widgets
from pyanc350.v2 import Positioner
import lakeshore
from scipy.optimize import curve_fit

def fitting_func(x, a, b, c, phi2, phi4):
    return a + b*np.sin((2*(x+phi2))/180*np.pi) + c*np.sin((4*(x+phi4))/180*np.pi)


def add_unique_postfix(fn):
    path, name = os.path.split(fn)
    name, ext = os.path.splitext(name)

    make_fn = lambda i: os.path.join(path, '%s_%04d%s' % (name, i, ext))

    for i in range(1, 1000):
        uni_fn = make_fn(i)
        if not os.path.exists(uni_fn):
            return uni_fn

    return None


def wait_time():
    # get filter & determine waittime
    filterN = daq.getInt('/%s/demods/0/order' % device)
    time_constant = daq.getDouble('/%s/demods/0/timeconstant' % device)
    if filterN == 1:
        waittime = 3 * time_constant
    elif filterN == 2:
        waittime = 4.7 * time_constant
    elif filterN == 3:
        waittime = 6.3 * time_constant
    elif filterN == 4:
        waittime = 7.8 * time_constant
    elif filterN == 5:
        waittime = 9.2 * time_constant
    elif filterN == 6:
        waittime = 11 * time_constant
    elif filterN == 7:
        waittime = 12 * time_constant
    elif filterN == 8:
        waittime = 13 * time_constant

    return waittime

f_conf = open(os.path.dirname(__file__)+ r'\..\..\Configuration.txt', "r")
conf_info = f_conf.read()
conf_info_split = conf_info.split('\n')
device_id = conf_info_split[0].split('\t')[1]
port_id = conf_info_split[1].split('\t')[1]
f_conf.close()

channel_name = ['/%s/demods/0/sample','/%s/demods/1/sample','/%s/demods/2/sample','/%s/demods/3/sample']

#device_id = 'DEV3425'
apilevel = 6
(daq, device, props) = ziutils.create_api_session(device_id, apilevel)

controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)
axis2_rot = newport.NewportESP301Axis(controller,2-1)
axis2_rot.enable()

axis1_rot = newport.NewportESP301Axis(controller,1-1)
axis1_rot.enable()

def find_sig_figs(number):
    dec = float(number)
    dec_to_str = str(dec)
    reverse = dec_to_str[::-1]
    return reverse.find('.')

def get_temp_values(start_temp, end_temp, step_size=0.):
    '''Returns an array of temperature values for a given starting temperature (starttemp), ending temperature (endtemp),
    and step size between temperature points (stepsize).
    starttemp: value to start list of temperatures at (in units of K)
    endtemp: maximum value for temperature to ramp up to (in units of K)
    stepsize: amount to increment temp by during ramp (in units of K)'''

    starttemp = float(start_temp)
    endtemp = float(end_temp)
    stepsize = float(step_size)
    temprange = endtemp-starttemp
    if stepsize == 0:
        return np.asarray(starttemp)
    else:
        roundsteps = np.rint([temprange/stepsize])[0]
        roundto = find_sig_figs(stepsize)
    if temprange/stepsize != roundsteps:
        newendtemp = starttemp + roundsteps*stepsize
        print('Warning: Specified end temperature value cannot be reached given the specified start temperature and step size. '
              +'End temperature will be rounded to the nearest step. New end temperature: '+str(newendtemp)+'K.')
        return np.around(np.linspace(starttemp, newendtemp, int(roundsteps)+1), roundto)
    else:
        return np.around(np.linspace(starttemp, endtemp, int(roundsteps)+1), roundto)

def get_angle_range(startangle, endangle, stepsize=0.):
    '''Returns an array of angle values for a given starting angle (startangle), ending angle (endangle),
    and step size between angle points (stepsize).
    startangle: value to start list of angles at (in units of deg)
    endangle: maximum value for angle to ramp up to (in units of deg)
    stepsize: amount to increment angle by during ramp (in units of deg)'''

    anglerange = endangle-startangle
    if stepsize == 0:
        return np.asarray(startangle)
    else:
        roundsteps = np.rint([anglerange/stepsize])[0]
    if anglerange/stepsize != roundsteps:
        newendangle = startangle + roundsteps*stepsize
        print('Warning: Specified end angle value cannot be reached given the specified start angle and step size. '
              +'End angle will be rounded to the nearest step. New end angle: '+str(newendangle)+' deg.')
        return np.linspace(startangle, newendangle, int(roundsteps)+1)
    else:
        return np.linspace(startangle, endangle, int(roundsteps)+1)

def get_z_values(start_z, end_z, step_size=1):
    '''Returns an array of temperature values for a given starting temperature (starttemp), ending temperature (endtemp),
    and step size between temperature points (stepsize).
    starttemp: value to start list of temperatures at (in units of K)
    endtemp: maximum value for temperature to ramp up to (in units of K)
    stepsize: amount to increment temp by during ramp (in units of K)'''

    startz=float(start_z)
    endz=float(end_z)
    stepsize=float(step_size)
    z_range = endz-startz
    if step_size == 0:
        return np.asarray(startz)
    else:
        roundsteps = np.rint([z_range/stepsize])[0]
    if z_range/stepsize != roundsteps:
        newendz = startz + roundsteps*stepsize
        print('Warning: Specified end z position cannot be reached given the specified start z and step size. '
              +'End z will be rounded to the nearest step. New end z: '+str(newend_z)+'um.')
        return np.linspace(startz, newendz, int(roundsteps)+1)
    else:
        return np.linspace(startz, endz, int(roundsteps)+1)

def get_input():
    global button
    while True:
        newValue = input('Stop? (Y/N) ')
        time.sleep(1)
        if (newValue == 'Y'):
            button = False
        if (button == False):
            break

def log_data_generator(filename_head):
    sf = 'S'
    while (sf!='Y') and (sf!='N'):
        sf = input('Save file? (Y/N)')
    if (sf == 'Y'):
        savefile = True
    elif (sf == 'N'):
        savefile = False
    if (savefile):
        sample_name = input('Sample? ')
        measurement_name = input('Measurement? ')
        time_stamp = str(datetime.datetime.now())
        print(time_stamp)
        time_stamp0 = time_stamp.replace(':','_')
        ax = {'x':0,'y':1,'z':2}
        anc = Positioner()
        print('x:',anc.getPosition(ax['x'])/1000)
        print('y:',anc.getPosition(ax['y'])/1000)
        print('z:',anc.getPosition(ax['z'])/1000)
        temp_stamp = input('Temperature (K)? ')
        field_stamp = input('Field (Oe)? ')
        objective_stamp = input('Objective? ')
        pump_wavelength = input('Pump wavelength (nm)? ')
        pump_power = input('Pump power (mW)? ')
        probe_wavelength = input('Probe wavelength (nm)? ')
        probe_power = input('Probe power (mW)? ')
        comments = input('Comments? ')
        filename = filename_head+ '\\'+str(sample_name) + '_'+ str(measurement_name) + '_' + time_stamp0+'.dat'
        file = open(filename,'a')
        file.write('Time'+'\t'+time_stamp+'\n')
        file.write('x_position'+'\t'+str(anc.getPosition(ax['x'])/1000)+'\n')
        file.write('y_position'+'\t'+str(anc.getPosition(ax['y'])/1000)+'\n')
        file.write('z_position'+'\t'+str(anc.getPosition(ax['z'])/1000)+'\n')
        file.write('Temperature'+'\t'+temp_stamp+'\t'+'K'+'\n')
        file.write('Field'+'\t'+field_stamp+'\t'+'Oe'+'\n')
        file.write('Objective'+'\t'+str(objective_stamp)+'\n')
        file.write('Pump wavelength'+'\t'+str(pump_wavelength)+'\t'+'nm'+'\n')
        file.write('Pump power'+'\t'+str(pump_power)+'\t'+'mW'+'\n')
        file.write('Probe wavelength'+'\t'+str(probe_wavelength)+'\t'+'nm'+'\n')
        file.write('Probe power'+'\t'+str(probe_power)+'\t'+'mW'+'\n')
        file.write('Comments'+'\t'+str(comments)+'\n')
        file.write('\n')
        file.close()
        anc.close()
        return [filename, savefile]
    else:
        return [0, savefile]


def log_data(filename_head, filename):
    ax = {'x': 0, 'y': 1, 'z': 2}
    anc = Positioner()
    time_stamp = str(datetime.datetime.now())
    filename1 = filename_head +'\\' +filename+'.dat'
    filename2 = add_unique_postfix(filename1)
    file = open(filename2,'a')
    file.write('Time'+'\t'+time_stamp+'\n')
    file.write('x_position'+'\t'+str(anc.getPosition(ax['x'])/1000)+'\n')
    file.write('y_position'+'\t'+str(anc.getPosition(ax['y'])/1000)+'\n')
    file.write('z_position'+'\t'+str(anc.getPosition(ax['z'])/1000)+'\n')
    file.write('Temperature'+'\t'+str(read_temperature())+'\t'+'K'+'\n')
    file.write('\n')
    file.close()
    anc.close()
    return filename2


def log_dataNoPos(filename_head, filename):

    time_stamp = str(datetime.datetime.now())
    filename1 = filename_head +'\\' +filename+'.dat'
    filename2 = add_unique_postfix(filename1)
    file = open(filename2,'a')
    file.write('Time'+'\t'+time_stamp+'\n')
    file.write('Temperature'+'\t'+str(read_temperature())+'\t'+'K'+'\n')
    file.write('\n')
    file.close()

    return filename2


def Lockin_time_record(channel_index, recording_time, time_constant, filename_head):
    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    logbook = log_data_generator(filename_head)
    filename = logbook[0]
    savefile = logbook[1]
    if (savefile):
        file = open(filename,'a')
        file.write("Time (s)"+'\t'+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
        file.close()

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

    time_record = np.array([])
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
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Demod x')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Demod y')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('R')
    draw_x, = ax1.plot([],'-o')
    draw_y, = ax2.plot([],'-o')
    draw_r, = ax3.plot([],'-o')
    fig.canvas.draw()
    fig.show()

    t_delay = 0
    tic = time.perf_counter()
    while (t_delay<recording_time):
        if button == False:
            break
        time.sleep(time_constant*4)
        toc = time.perf_counter()
        t_delay = toc - tic
        sample = daq.getSample(channel_name[channel_index-1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        time_record = np.append(time_record, t_delay)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        #Plot
        draw_x.set_data(time_record-time_record[0],demod_x)
        draw_y.set_data(time_record-time_record[0],demod_y)
        draw_r.set_data(time_record-time_record[0],demod_r)
        ax1.relim()
        ax1.autoscale()
        ax2.relim()
        ax2.autoscale()
        ax3.relim()
        ax3.autoscale()
        fig.canvas.draw()
        fig.canvas.flush_events()

        #Save file
        if (savefile):
            file = open(filename,'a')
            file.write(format(t_delay-time_record[0], '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
            file.close()
    thread1.join()


def Find_balance_angle(incident_pol_angle, axis_index, start_pos, step_size, num_of_steps, go_back, channel_index, time_constant, filename_head, filename):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    filename = log_data(filename_head, filename)
    savefile =1
    print('Balance for', incident_pol_angle, 'incident polarization')
    if (savefile):
        file = open(filename,'a')
        file.write('Balance for '+str(incident_pol_angle)+' incident polarization'+'\n')
        file.write('\n')
        file.write("Angle (deg)"+'\t'+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
        file.close()

    #Scan parameter
    end_pos = start_pos + step_size * (num_of_steps-1)
    scan_range = np.arange(incident_pol_angle+start_pos, incident_pol_angle+start_pos + step_size * num_of_steps, step_size)
    axis_rot = newport.NewportESP301Axis(controller,axis_index-1)
    axis_rot.enable()
    go_to_balance_angle = input('Go to balance angle? (True/False) ')
    real_start = start_pos - go_back

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

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

    #Scan
    for pos in scan_range:
        if button == False:
            break
        if (pos == start_pos):
            axis_rot.move(pos-go_back,absolute=True)
            while (axis_rot.is_motion_done==False):
                pass
        axis_rot.move(pos,absolute=True)
        while (axis_rot.is_motion_done==False):
            pass
        time.sleep(time_constant*4)
        sample = daq.getSample(channel_name[int(channel_index-1)] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        position = np.append(position, pos)
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

        #Save file
        if (savefile):
            file = open(filename,'a')
            file.write(format(pos, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
            file.close()

    if (button == True):
        k = np.polyfit(scan_range[1:], demod_x[1:], 1)
        print(k[0])
        balance_angle = -k[1]/k[0]
        pos_ref = np.linspace(incident_pol_angle+start_pos, incident_pol_angle+end_pos, 1000)

        ax1.plot(pos_ref, k[0]*pos_ref+k[1],'C1')
        fig.canvas.draw()
        fig.canvas.flush_events()

        if (savefile):
            file = open(filename,'a')
            file.write('Balance angle = '+str(balance_angle)+' deg'+'\n')
            file.close()

        print('Balance angle = '+str(balance_angle)+' deg')
        print('Current angle =')


        if (go_to_balance_angle):
            dh = display(str(axis_rot.position), display_id=True)
            #axis_rot.move(balance_angle-go_back,absolute=True)
            axis_rot.move(real_start,absolute=True)
            while (axis_rot.is_motion_done==False):
                    dh.update(str(axis_rot.position))
            axis_rot.move(balance_angle,absolute=True)
            while (axis_rot.is_motion_done==False):
                    dh.update(str(axis_rot.position))
    thread1.join()

def Rotate_one_HWP(axis_index, start_pos, step_size, num_of_steps, go_back, channel_index,filename_head, filename):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    filename = log_data(filename_head, filename)
    savefile = 1
    if (savefile):
        file = open(filename, 'a')
        file.write("Angle_1 (deg)" + '\t' + "Angle_2 (deg)" + '\t' + "Demod x" + '\t' + "Demod y" + '\t' + "R" + '\n')
        file.close()
    waittime = wait_time()
    end_pos = start_pos + step_size * (num_of_steps - 1)
    scan_range = np.arange(start_pos, start_pos + step_size* num_of_steps, step_size)
    axis_index_1=1
    axis_index_2=2
    axis_rot_1 = newport.NewportESP301Axis(controller, axis_index_1 - 1)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller, axis_index_2 - 1)
    axis_rot_2.enable()

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

    # Quantities to be measured
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    fig = plt.figure(figsize=(8, 10))
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
    draw_x, = ax1.plot([], '-o')
    draw_y, = ax2.plot([], '-o')
    draw_r, = ax3.plot([], '-o')
    fig.canvas.draw()
    fig.show()

    # Scan
    for i in np.arange(0, len(scan_range), 1):
        if button == False:
            break
        pos = scan_range[i]
        if (pos == start_pos):
            if axis_index==1:
                axis_rot_1.move(pos - go_back, absolute=True)
            elif axis_index==2:
                axis_rot_2.move(pos - go_back, absolute=True)
            while (axis_rot_1.is_motion_done == False) or (axis_rot_2.is_motion_done == False):
                pass
        if axis_index==1:
            axis_rot_1.move(pos, absolute=True)
        elif axis_index==2:
            axis_rot_2.move(pos, absolute=True)
        while (axis_rot_1.is_motion_done == False) or (axis_rot_2.is_motion_done == False):
            pass
        time.sleep(waittime)
        sample = daq.getSample(channel_name[channel_index - 1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        if axis_index==1:
            position = np.append(position, 2 * axis_rot_1.position)
        elif axis_index==2:
            position = np.append(position, 2 * axis_rot_2.position)

        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # Plot
        draw_x.set_data(position, demod_x)
        draw_y.set_data(position, demod_y)
        draw_r.set_data(position, demod_r)
        ax1.relim()
        ax1.autoscale()
        ax2.relim()
        ax2.autoscale()
        ax3.relim()
        ax3.autoscale()
        fig.canvas.draw()
        fig.canvas.flush_events()

        # Save file
        if (savefile):
            file = open(filename, 'a')
            file.write(
                format(axis_rot_1.position, '.15f') + "\t" + format(axis_rot_2.position, '.15f') + "\t" + format(x,
                                                                                                                 '.15f') + '\t' + format(
                    y, '.15f') + '\t' + format(r, '.15f') + '\n')
            file.close()

    popt, popv = curve_fit(fitting_func, position, demod_x, p0=[0, 1e-5,0, 0, 0 ], bounds=([-0.1, 0, 0, 0, 0], [0.1, 0.1, 1e-5, 180, 90]))
    print(popt)
    pos_ref = np.linspace(2 * start_pos, 2 * end_pos, 1000)

    ax1.plot(pos_ref, fitting_func(pos_ref, popt[0], popt[1], popt[2], popt[3], popt[4]), 'C1')
    fig.canvas.draw()
    fig.canvas.flush_events()
    thread1.join()

def Corotate_measurement(num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, step_size_2, go_back_2, channel_index, time_constant, filename_head, filename, showplot):
    # ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    filename = log_data(filename_head, filename)
    savefile = 1
    if (savefile):
        file = open(filename, 'a')
        file.write("Time (s)"+"\t"+"Angle_1 (deg)" + '\t' + "Angle_2 (deg)" + '\t' + "Demod x" + '\t' + "Demod y" + '\t' + "R"+'\t' + "Demod x" + '\t' + "Demod y" + '\t' + "R" +'\n')
        file.close()
    waittime = wait_time()
    end_pos_1 = start_pos_1 + step_size_1 * (num_of_steps - 1)
    end_pos_2 = start_pos_2 + step_size_2 * (num_of_steps - 1)
    scan_range_1 = np.arange(start_pos_1, start_pos_1 + step_size_1 * num_of_steps, step_size_1)
    scan_range_2 = np.arange(start_pos_2, start_pos_2 + step_size_2 * num_of_steps, step_size_2)
    axis_rot_1 = newport.NewportESP301Axis(controller, axis_index_1 - 1)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller, axis_index_2 - 1)
    axis_rot_2.enable()
    start_pos = start_pos_1
    end_pos = start_pos_1 + step_size_1 * (num_of_steps - 1)

    #global button
    #button = True
    #thread1 = threading.Thread(target=get_input)
    #thread1.start()

    # Quantities to be measured
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    if showplot:
        fig = plt.figure(figsize=(8, 10))
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
        draw_x, = ax1.plot([], '-o')
        draw_y, = ax2.plot([], '-o')
        draw_r, = ax3.plot([], '-o')
        fig.canvas.draw()
        fig.show()

    t0=time.time();


    # Scan
    for i in np.arange(0, len(scan_range_1), 1):
        #if button == False:
    #        break
        pos_1 = scan_range_1[i]
        pos_2 = scan_range_2[i]
        if (pos_1 == start_pos_1):
            axis_rot_1.move(pos_1 - go_back_1, absolute=True)
            axis_rot_2.move(pos_2 - go_back_2, absolute=True)
            while True:
                time.sleep(0.03)
                try:
                    if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                        break
                except ValueError:
                    pass
        axis_rot_1.move(pos_1, absolute=True)
        time.sleep(0.03)
        axis_rot_2.move(pos_2, absolute=True)
        while True:
            time.sleep(0.03)
            try:
                if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                    break
            except ValueError:
                pass
        time.sleep(waittime)
        sample = daq.getSample(channel_name[channel_index - 1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]

        sample1 = daq.getSample(channel_name[3 - 1] % device)
        sample1["R"] = np.abs(sample1["x"] + 1j * sample1["y"])
        x1 = sample1["x"][0]
        y1= sample1["y"][0]
        r1 = sample1["R"][0]



        position = np.append(position, 2 * axis_rot_1.position)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # Plot
        if showplot:
            draw_x.set_data(position, demod_x)
            draw_y.set_data(position, demod_y)
            draw_r.set_data(position, demod_r)
            ax1.relim()
            ax1.autoscale()
            ax2.relim()
            ax2.autoscale()
            ax3.relim()
            ax3.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

        # Save file
        t=time.time()-t0
        if (savefile):
            file = open(filename, 'a')
            file.write(
                str(t)+ "\t" +format(axis_rot_1.position, '.15f')[0:4] + "\t" + format(axis_rot_2.position, '.15f')[0:4]  + "\t" + format(x, '.15f') + '\t' + format(y, '.15f') + '\t' + format(r, '.15f') + "\t" + format(x1, '.15f') + '\t' + format(y1, '.15f') + '\t' + format(r1, '.15f')+'\n')
            file.close()

def Corotate_measurementNoPos(num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, step_size_2, go_back_2, channel_index, time_constant, filename_head, filename, showplot):
    # ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)
    daq.setDouble('/%s/demods/0/timeconstant' % device, time_constant)
    filename = log_dataNoPos(filename_head, filename)
    savefile = 1
    if (savefile):
        file = open(filename, 'a')
        file.write("Time (s)"+"\t"+"Angle_1 (deg)" + '\t' + "Angle_2 (deg)" + '\t' + "Demod x" + '\t' + "Demod y" + '\t' + "R"+'\t' + "Demod x" + '\t' + "Demod y" + '\t' + "R" +'\n')
        file.close()
    waittime = wait_time()
    end_pos_1 = start_pos_1 + step_size_1 * (num_of_steps - 1)
    end_pos_2 = start_pos_2 + step_size_2 * (num_of_steps - 1)
    scan_range_1 = np.arange(start_pos_1, start_pos_1 + step_size_1 * num_of_steps, step_size_1)
    scan_range_2 = np.arange(start_pos_2, start_pos_2 + step_size_2 * num_of_steps, step_size_2)
    axis_rot_1 = newport.NewportESP301Axis(controller, axis_index_1 - 1)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller, axis_index_2 - 1)
    axis_rot_2.enable()
    start_pos = start_pos_1
    end_pos = start_pos_1 + step_size_1 * (num_of_steps - 1)

    #global button
    #button = True
    #thread1 = threading.Thread(target=get_input)
    #thread1.start()

    # Quantities to be measured
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    if showplot:
        fig = plt.figure(figsize=(8, 10))
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
        draw_x, = ax1.plot([], '-o')
        draw_y, = ax2.plot([], '-o')
        draw_r, = ax3.plot([], '-o')
        fig.canvas.draw()
        fig.show()

    t0=time.time();


    # Scan
    for i in np.arange(0, len(scan_range_1), 1):
        #if button == False:
    #        break
        pos_1 = scan_range_1[i]
        pos_2 = scan_range_2[i]
        if (pos_1 == start_pos_1):
            axis_rot_1.move(pos_1 - go_back_1, absolute=True)
            axis_rot_2.move(pos_2 - go_back_2, absolute=True)
            while True:
                time.sleep(0.03)
                try:
                    if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                        break
                except ValueError:
                    pass
        axis_rot_1.move(pos_1, absolute=True)
        time.sleep(0.03)
        axis_rot_2.move(pos_2, absolute=True)
        while True:
            time.sleep(0.03)
            try:
                if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                    break
            except ValueError:
                pass
        time.sleep(waittime)
        sample = daq.getSample(channel_name[channel_index - 1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]

        sample1 = daq.getSample(channel_name[3 - 1] % device)
        sample1["R"] = np.abs(sample1["x"] + 1j * sample1["y"])
        x1 = sample1["x"][0]
        y1= sample1["y"][0]
        r1 = sample1["R"][0]



        position = np.append(position, 2 * axis_rot_1.position)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # Plot
        if showplot:
            draw_x.set_data(position, demod_x)
            draw_y.set_data(position, demod_y)
            draw_r.set_data(position, demod_r)
            ax1.relim()
            ax1.autoscale()
            ax2.relim()
            ax2.autoscale()
            ax3.relim()
            ax3.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

        # Save file
        t=time.time()-t0
        if (savefile):
            file = open(filename, 'a')
            file.write(
                str(t)+ "\t" +format(axis_rot_1.position, '.15f')[0:4] + "\t" + format(axis_rot_2.position, '.15f')[0:4]  + "\t" + format(x, '.15f') + '\t' + format(y, '.15f') + '\t' + format(r, '.15f') + "\t" + format(x1, '.15f') + '\t' + format(y1, '.15f') + '\t' + format(r1, '.15f')+'\n')
            file.close()

def Corotate_measurement_fit(num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, step_size_2, go_back_2, channel_index, time_constant, filename_head, filename):
    # ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    filename = log_data(filename_head, filename)
    savefile = 1
    if (savefile):
        file = open(filename, 'a')
        file.write("Time (s)"+"\t"+"Angle_1 (deg)" + '\t' + "Angle_2 (deg)" + '\t' + "Demod x" + '\t' + "Demod y" + '\t' + "R" + '\n')
        file.close()
    waittime = wait_time()
    end_pos_1 = start_pos_1 + step_size_1 * (num_of_steps - 1)
    end_pos_2 = start_pos_2 + step_size_2 * (num_of_steps - 1)
    scan_range_1 = np.arange(start_pos_1, start_pos_1 + step_size_1 * num_of_steps, step_size_1)
    scan_range_2 = np.arange(start_pos_2, start_pos_2 + step_size_2 * num_of_steps, step_size_2)
    axis_rot_1 = newport.NewportESP301Axis(controller, axis_index_1 - 1)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller, axis_index_2 - 1)
    axis_rot_2.enable()
    start_pos = start_pos_1
    end_pos = start_pos_1 + step_size_1 * (num_of_steps - 1)

    #global button
    #button = True
    #thread1 = threading.Thread(target=get_input)
    #thread1.start()

    # Quantities to be measured
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    fig = plt.figure(figsize=(8, 10))
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
    draw_x, = ax1.plot([], '-o')
    draw_y, = ax2.plot([], '-o')
    draw_r, = ax3.plot([], '-o')
    fig.canvas.draw()
    fig.show()

    t0=time.time();


    # Scan
    for i in np.arange(0, len(scan_range_1), 1):
        #if button == False:
    #        break
        pos_1 = scan_range_1[i]
        pos_2 = scan_range_2[i]
        if (pos_1 == start_pos_1):
            axis_rot_1.move(pos_1 - go_back_1, absolute=True)
            axis_rot_2.move(pos_2 - go_back_2, absolute=True)
            while True:
                time.sleep(0.03)
                try:
                    if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                        break
                except ValueError:
                    pass
        axis_rot_1.move(pos_1, absolute=True)
        time.sleep(0.03)
        axis_rot_2.move(pos_2, absolute=True)
        while True:
            time.sleep(0.03)
            try:
                if axis_rot_1.is_motion_done==True and axis_rot_2.is_motion_done==True:
                    break
            except ValueError:
                pass
        time.sleep(waittime)
        sample = daq.getSample(channel_name[channel_index - 1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        position = np.append(position, 2 * axis_rot_1.position)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        # Plot
        draw_x.set_data(position, demod_x)
        draw_y.set_data(position, demod_y)
        draw_r.set_data(position, demod_r)
        ax1.relim()
        ax1.autoscale()
        ax2.relim()
        ax2.autoscale()
        ax3.relim()
        ax3.autoscale()
        fig.canvas.draw()
        fig.canvas.flush_events()

        # Save file
        t=time.time()-t0
        if (savefile):
            file = open(filename, 'a')
            file.write(
                str(t)+ "\t" +format(axis_rot_1.position, '.15f')[0:4] + "\t" + format(axis_rot_2.position, '.15f')[0:4]  + "\t" + format(x,
                                                                                                                 '.15f') + '\t' + format(
                    y, '.15f') + '\t' + format(r, '.15f') + '\n')
            file.close()

    popt, popv = curve_fit(fitting_func, position, demod_x,  bounds=([-0.1, 0, 0, 0, 0], [0.1, 0.1, 1e-3, 180, 90]))

    print(popt)
    pos_ref = np.linspace(2 * start_pos, 2 * end_pos, 1000)

    ax1.plot(pos_ref, fitting_func(pos_ref, popt[0], popt[1], popt[2], popt[3], popt[4]), 'C1')
    fig.canvas.draw()
    fig.canvas.flush_events()
    #thread1.join()


def Pump_probe(axis_index, start_pos, step_size, num_of_steps, go_back, channel_index, time_constant, filename_head):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)


    #Filename & title
    logbook = log_data_generator(filename_head)
    filename = logbook[0]
    savefile = logbook[1]
    if (savefile):
        file = open(filename,'a')
        file.write("Position (mm)"+"\t"+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
        file.close()

    #Scan parameter
    scan_range = np.arange(start_pos, start_pos + step_size * num_of_steps, step_size)
    axis_delay = newport.NewportESP301Axis(controller,axis_index-1)
    axis_delay.enable()

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

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
    ax1.set_xlabel('Position (mm)')
    ax1.set_ylabel('Demod x')
    ax2.set_xlabel('Position (mm)')
    ax2.set_ylabel('Demod y')
    ax3.set_xlabel('Position (mm)')
    ax3.set_ylabel('R')
    draw_x, = ax1.plot([],'-o')
    draw_y, = ax2.plot([],'-o')
    draw_r, = ax3.plot([],'-o')
    fig.canvas.draw()
    fig.show()

    #Scan
    for pos in scan_range:
        if button == False:
            break
        if (pos == start_pos):
            axis_delay.move(pos-go_back,absolute=True)
            while (axis_delay.is_motion_done==False):
                pass
        axis_delay.move(pos,absolute=True)
        while (axis_delay.is_motion_done==False):
            pass
        time.sleep(time_constant*4)
        sample = daq.getSample('/%s/demods/0/sample' % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        position = np.append(position, pos)
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

        #Save file
        if (savefile):
            file = open(filename,'a')
            file.write(format(pos, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
            file.close()
    thread1.join()


def MappingY(filename_head):
    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    #Filename & title
    logbook = log_data_generator(filename_head)
    filename = logbook[0]
    savefile = logbook[1]
    if (savefile):
        file = open(filename,'a')
        file.write("x coordinate (um)"+"\t"+"y coordinate (um)"+"\t"+"z coordinate (um)"+"\t"+"Demod x"+"\t"+"Demod y"+"\t"+"R"+"\n")
        file.close()

    #Attocube initialization
    ax = {'x':0,'y':1,'z':2}
    anc = Positioner()

    #Scan parameter
    x_start = float(input('x_start_position? '))
    x_step = float(input('x_step? '))
    x_num = int(input('x_numbers? '))
    x_end = x_start + x_step * x_num

    y_start = float(input('y_start_position? '))
    y_step = float(input('y_step? '))
    y_num = int(input('y_numbers? '))
    y_end = y_start + y_step * y_num

    x_range = np.arange(x_start, x_end, x_step)
    y_range = np.arange(y_start, y_end, y_step)
    x_plot_range = np.arange(x_start, x_end + x_step, x_step)
    y_plot_range = np.arange(y_start, y_end + y_step, y_step)
    x_tor = float(input('x_tolerance? '))
    y_tor = float(input('y_tolerance? '))
    time_constant = float(input('time_constant? '))   #unit: second
    go_back = float(input('go_back? '))          #preventing hysteresis
    channel_index = int(input('channel_index? '))
    R_channel_index = int(input('R_channel_index? '))

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

    x_position = np.array([])
    y_position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])
    demod_x_R = np.array([])
    demod_y_R = np.array([])
    demod_r_R = np.array([])

    fig = plt.figure(figsize=(5,15))
    gs = fig.add_gridspec(4, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[3, 0])
    ax1.set_xlabel('x (um)')
    ax1.set_ylabel('y (um)')
    ax1.set_xticks(x_plot_range)
    ax1.set_xticklabels(x_plot_range,rotation=90)
    ax1.set_yticks(y_plot_range)
    ax2.set_xlabel('x (um)')
    ax2.set_ylabel('y (um)')
    ax2.set_xticks(x_plot_range)
    ax2.set_xticklabels(x_plot_range,rotation=90)
    ax2.set_yticks(y_plot_range)
    ax3.set_xlabel('x (um)')
    ax3.set_ylabel('y (um)')
    ax3.set_xticks(x_plot_range)
    ax3.set_xticklabels(x_plot_range,rotation=90)
    ax3.set_yticks(y_plot_range)
    ax4.set_xlabel('x (um)')
    ax4.set_ylabel('y (um)')
    ax4.set_xticks(x_plot_range)
    ax4.set_xticklabels(x_plot_range,rotation=90)
    ax4.set_yticks(y_plot_range)

    demod_x0 = np.zeros((y_num, x_num))
    demod_y0 = np.zeros((y_num, x_num))
    demod_r0 = np.zeros((y_num, x_num))
    demod_r_R0 = np.zeros((y_num, x_num))
    x_coordinates = np.arange(x_start, x_end, x_step)
    y_coordinates = np.arange(y_start, y_end, y_step)
    X_coor, Y_coor = np.meshgrid(x_coordinates, y_coordinates)
    extent=[x_start, x_end, y_start, y_end]

    pos1 = ax1.imshow(demod_x0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos1, ax=ax1)
    pos2 = ax2.imshow(demod_y0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos2, ax=ax2)
    pos3 = ax3.imshow(demod_r0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos3, ax=ax3)
    pos4 = ax4.imshow(demod_r_R0,cmap="bwr",origin='lower',extent=extent)
    fig.colorbar(pos4, ax=ax4)
    fig.canvas.draw()
    fig.tight_layout()
    fig.show()



    for y_pos in y_range:
        if button == False:
            break
        for x_pos in x_range:
            if button == False:
                break
            if (x_pos == x_start) and (y_pos == y_start):
                x_target = x_pos-go_back
                y_target = y_pos-go_back
                anc.moveAbsolute(ax['x'], int(x_target*1000))
                anc.moveAbsolute(ax['y'], int(y_target*1000))
                x_error = np.abs(x_target-anc.getPosition(ax['x'])/1000)
                y_error = np.abs(y_target-anc.getPosition(ax['y'])/1000)
                while (x_error >= x_tor) or (y_error >= y_tor):
                    if button == False:
                        break
                    time.sleep(1)
                    #clear_output(wait=True)
                    x_error = np.abs(x_target-anc.getPosition(ax['x'])/1000)
                    y_error = np.abs(y_target-anc.getPosition(ax['y'])/1000)
                    #print(x_error)
                    #print(y_error)
                    #print(anc.getStatus(ax['x']))
                    #print(anc.getStatus(ax['y']))
                    if (x_error >= x_tor):
                        anc.moveAbsolute(ax['x'], int(x_target*1000))
                    if (y_error >= y_tor):
                        anc.moveAbsolute(ax['y'], int(y_target*1000))
                if button == False:
                    break
            anc.moveAbsolute(ax['x'], int(x_pos*1000))
            x_error = np.abs(x_pos-anc.getPosition(ax['x'])/1000)
            while (x_error >= x_tor):
                if button == False:
                    break
                time.sleep(1)
                #clear_output(wait=True)
                x_error = np.abs(x_pos-anc.getPosition(ax['x'])/1000)
                #print(x_error)
                if (x_error >= x_tor):
                    anc.moveAbsolute(ax['x'], int(x_pos*1000))
            if button == False:
                break
            anc.moveAbsolute(ax['y'], int(y_pos*1000))
            y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
            while (y_error >= y_tor):
                if button == False:
                    break
                time.sleep(1)
                #clear_output(wait=True)
                y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
                #print(y_error)
                if (y_error >= y_tor):
                    anc.moveAbsolute(ax['y'], int(y_pos*1000))
            if button == False:
                break
            print("moved to "+str(anc.getPosition(ax['x'])/1000)+","+str(anc.getPosition(ax['y'])/1000))
            time.sleep(time_constant*4)
            x_position = np.append(x_position, anc.getPosition(ax['x'])/1000)
            y_position = np.append(y_position, anc.getPosition(ax['y'])/1000)
            sample = daq.getSample(channel_name[channel_index-1] % device)
            sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
            x = sample["x"][0]
            y = sample["y"][0]
            r = sample["R"][0]
            demod_x = np.append(demod_x, x)
            demod_y = np.append(demod_y, y)
            demod_r = np.append(demod_r, r)

            sample_R = daq.getSample(channel_name[R_channel_index-1] % device)
            sample_R["R"] = np.abs(sample_R["x"] + 1j * sample_R["y"])
            x_R = sample_R["x"][0]
            y_R = sample_R["y"][0]
            r_R = sample_R["R"][0]
            demod_x_R = np.append(demod_x_R, x_R)
            demod_y_R = np.append(demod_y_R, y_R)
            demod_r_R = np.append(demod_r_R, r_R)

             #Plot
            length = len(demod_x)
            y_num0 = length//x_num
            x_num0 = length-y_num0*x_num
            demod_x0[:y_num0, :] = np.reshape(demod_x[:y_num0*x_num], (y_num0, x_num))
            demod_y0[:y_num0, :] = np.reshape(demod_y[:y_num0*x_num], (y_num0, x_num))
            demod_r0[:y_num0, :] = np.reshape(demod_r[:y_num0*x_num], (y_num0, x_num))
            demod_r_R0[:y_num0, :] = np.reshape(demod_r_R[:y_num0*x_num], (y_num0, x_num))
            if (y_num0 < y_num):
                demod_x0[y_num0, :x_num0] = demod_x[y_num0*x_num:length]
                demod_y0[y_num0, :x_num0] = demod_y[y_num0*x_num:length]
                demod_r0[y_num0, :x_num0] = demod_r[y_num0*x_num:length]
                demod_r_R0[y_num0, :x_num0] = demod_r_R[y_num0*x_num:length]

            pos1.set_data(demod_x0)
            pos1.set_clim(vmin = demod_x0.min(), vmax = demod_x0.max())
            pos2.set_data(demod_y0)
            pos2.set_clim(vmin = demod_y0.min(), vmax = demod_y0.max())
            pos3.set_data(demod_r0)
            pos3.set_clim(vmin = demod_r0.min(), vmax = demod_r0.max())
            pos4.set_data(demod_r_R0)
            pos4.set_clim(vmin = demod_r_R0.min(), vmax = demod_r_R0.max())
            fig.canvas.draw()
            fig.canvas.flush_events()


            if (savefile):
                file = open(filename,'a')
                file.write(format(x_position[len(x_position)-1], '.15f')+"\t"+format(y_position[len(y_position)-1], '.15f')+"\t"+format(anc.getPosition(ax['z'])/1000, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\t'+format(r_R, '.15f')+'\n')
                file.close()

    anc.close()
    thread1.join()


def Mapping(filename_head,filename, logging, x_start, x_end, y_start, y_end, step, tol, channel_index, R_channel_index, time_constant ):
    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    # Filename & title

    if logging:
        logbook = log_data_generator(filename_head)
        filename = logbook[0]
        savefile = logbook[1]
    else:
        savefile = 1
        filename=log_data(filename_head, filename)

    if (savefile):
        file = open(filename, 'a')
        file.write(
            "time (s)" + "\t" +"x coordinate (um)" + "\t" + "y coordinate (um)" + "\t" + "z coordinate (um)" + "\t" + "Demod x" + "\t" + "Demod y" + "\t" + "R" + "\n")
        file.close()

    waittime = wait_time()

    # Attocube initialization
    ax = {'x': 0, 'y': 1, 'z': 2}
    anc = Positioner()

    # Scan parameter

    x_step=step
    y_step=step

    x_num = int( abs((x_end-x_start)/x_step+1))
    y_num = int(abs((y_end-y_start)/y_step+1))


    x_range = np.arange(x_start, x_end+x_step, x_step)
    y_range = np.arange(y_start, y_end+x_step, y_step)
    x_plot_range = np.arange(x_start, x_end + x_step, x_step)
    y_plot_range = np.arange(y_start, y_end + y_step, y_step)
    x_tor =tol
    y_tor = tol
    go_back = 10  # preventing hysteresis

    global button
    button = True
    #thread1 = threading.Thread(target=get_input)
    #thread1.start()

    x_position = np.array([])
    y_position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])
    demod_x_R = np.array([])
    demod_y_R = np.array([])
    demod_r_R = np.array([])

    fig = plt.figure(figsize=(5, 15))
    gs = fig.add_gridspec(4, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[3, 0])
    ax1.set_xlabel('x (um)')
    ax1.set_ylabel('y (um)')
    ax1.set_xticks(x_plot_range)
    ax1.set_xticklabels(x_plot_range, rotation=90)
    ax1.set_yticks(y_plot_range)
    ax2.set_xlabel('x (um)')
    ax2.set_ylabel('y (um)')
    ax2.set_xticks(x_plot_range)
    ax2.set_xticklabels(x_plot_range, rotation=90)
    ax2.set_yticks(y_plot_range)
    ax3.set_xlabel('x (um)')
    ax3.set_ylabel('y (um)')
    ax3.set_xticks(x_plot_range)
    ax3.set_xticklabels(x_plot_range, rotation=90)
    ax3.set_yticks(y_plot_range)
    ax4.set_xlabel('x (um)')
    ax4.set_ylabel('y (um)')
    ax4.set_xticks(x_plot_range)
    ax4.set_xticklabels(x_plot_range, rotation=90)
    ax4.set_yticks(y_plot_range)

    demod_x0 = np.zeros((y_num, x_num))
    demod_y0 = np.zeros((y_num, x_num))
    demod_r0 = np.zeros((y_num, x_num))
    demod_r_R0 = np.zeros((y_num, x_num))
    x_coordinates = np.arange(x_start, x_end, x_step)
    y_coordinates = np.arange(y_start, y_end, y_step)
    X_coor, Y_coor = np.meshgrid(x_coordinates, y_coordinates)
    extent = [x_start, x_end, y_start, y_end]

    pos1 = ax1.imshow(demod_x0, cmap="bwr", origin='lower', extent=extent, norm=colors.TwoSlopeNorm(0))
    fig.colorbar(pos1, ax=ax1)
    pos2 = ax2.imshow(demod_y0, cmap="bwr", origin='lower', extent=extent, norm=colors.TwoSlopeNorm(0))
    fig.colorbar(pos2, ax=ax2)
    pos3 = ax3.imshow(demod_r0, cmap="bwr", origin='lower', extent=extent, norm=colors.TwoSlopeNorm(0))
    fig.colorbar(pos3, ax=ax3)
    pos4 = ax4.imshow(demod_r_R0, cmap="bwr", origin='lower', extent=extent)
    fig.colorbar(pos4, ax=ax4)
    fig.canvas.draw()
    fig.tight_layout()
    fig.show()
    t0 = time.time()
    for y_pos in y_range:
        if button == False:
            break
        for x_pos in x_range:
            if button == False:
                break
            if (x_pos == x_start) and (y_pos == y_start):
                x_target = x_pos - go_back
                y_target = y_pos - go_back
                anc.moveAbsolute(ax['x'], int(x_target * 1000))
                anc.moveAbsolute(ax['y'], int(y_target * 1000))
                x_error = np.abs(x_target - anc.getPosition(ax['x']) / 1000)
                y_error = np.abs(y_target - anc.getPosition(ax['y']) / 1000)
                while (x_error >= x_tor) or (y_error >= y_tor):
                    if button == False:
                        break
                    time.sleep(1)
                    # clear_output(wait=True)
                    x_error = np.abs(x_target - anc.getPosition(ax['x']) / 1000)
                    y_error = np.abs(y_target - anc.getPosition(ax['y']) / 1000)
                    # print(x_error)
                    # print(y_error)
                    # print(anc.getStatus(ax['x']))
                    # print(anc.getStatus(ax['y']))
                    if (x_error >= x_tor):
                        anc.moveAbsolute(ax['x'], int(x_target * 1000))
                    if (y_error >= y_tor):
                        anc.moveAbsolute(ax['y'], int(y_target * 1000))
                if button == False:
                    break
            anc.moveAbsolute(ax['x'], int(x_pos * 1000))
            x_error = np.abs(x_pos - anc.getPosition(ax['x']) / 1000)
            while (x_error >= x_tor):
                if button == False:
                    break
                time.sleep(1)
                # clear_output(wait=True)
                x_error = np.abs(x_pos - anc.getPosition(ax['x']) / 1000)
                # print(x_error)
                if (x_error >= x_tor):
                    anc.moveAbsolute(ax['x'], int(x_pos * 1000))
            if button == False:
                break
            anc.moveAbsolute(ax['y'], int(y_pos * 1000))
            y_error = np.abs(y_pos - anc.getPosition(ax['y']) / 1000)
            while (y_error >= y_tor):
                if button == False:
                    break
                time.sleep(1)
                # clear_output(wait=True)
                y_error = np.abs(y_pos - anc.getPosition(ax['y']) / 1000)
                # print(y_error)
                if (y_error >= y_tor):
                    anc.moveAbsolute(ax['y'], int(y_pos * 1000))
            if button == False:
                break
            print("moved to " + str(anc.getPosition(ax['x']) / 1000) + "," + str(anc.getPosition(ax['y']) / 1000))
            time.sleep(waittime)
            x_position = np.append(x_position, anc.getPosition(ax['x']) / 1000)
            y_position = np.append(y_position, anc.getPosition(ax['y']) / 1000)
            sample = daq.getSample(channel_name[channel_index - 1] % device)
            sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
            x = sample["x"][0]
            y = sample["y"][0]
            r = sample["R"][0]
            demod_x = np.append(demod_x, x)
            demod_y = np.append(demod_y, y)
            demod_r = np.append(demod_r, r)

            sample_R = daq.getSample(channel_name[R_channel_index - 1] % device)
            sample_R["R"] = np.abs(sample_R["x"] + 1j * sample_R["y"])
            x_R = sample_R["x"][0]
            y_R = sample_R["y"][0]
            r_R = sample_R["R"][0]
            demod_x_R = np.append(demod_x_R, x_R)
            demod_y_R = np.append(demod_y_R, y_R)
            demod_r_R = np.append(demod_r_R, r_R)

            # Plot
            length = len(demod_x)
            y_num0 = length // x_num
            x_num0 = length - y_num0 * x_num
            demod_x0[:y_num0, :] = np.reshape(demod_x[:y_num0 * x_num], (y_num0, x_num))
            demod_y0[:y_num0, :] = np.reshape(demod_y[:y_num0 * x_num], (y_num0, x_num))
            demod_r0[:y_num0, :] = np.reshape(demod_r[:y_num0 * x_num], (y_num0, x_num))
            demod_r_R0[:y_num0, :] = np.reshape(demod_r_R[:y_num0 * x_num], (y_num0, x_num))
            if (y_num0 < y_num):
                demod_x0[y_num0, :x_num0] = demod_x[y_num0 * x_num:length]
                demod_y0[y_num0, :x_num0] = demod_y[y_num0 * x_num:length]
                demod_r0[y_num0, :x_num0] = demod_r[y_num0 * x_num:length]
                demod_r_R0[y_num0, :x_num0] = demod_r_R[y_num0 * x_num:length]

            pos1.set_data(demod_x0)
            pos1.set_clim(vmin=demod_x0.min(), vmax=demod_x0.max())
            pos2.set_data(demod_y0)
            pos2.set_clim(vmin=demod_y0.min(), vmax=demod_y0.max())
            pos3.set_data(demod_r0)
            pos3.set_clim(vmin=demod_r0.min(), vmax=demod_r0.max())
            pos4.set_data(demod_r_R0)
            pos4.set_clim(vmin=demod_r_R0.min(), vmax=demod_r_R0.max())
            fig.canvas.draw()
            fig.canvas.flush_events()

            if (savefile):
                # time
                ts = time.time() - t0;
                file = open(filename, 'a')
                file.write(format(ts, '.15f') + "\t"+
                    format(x_position[len(x_position) - 1], '.15f') + "\t" + format(y_position[len(y_position) - 1],
                                                                                    '.15f') + "\t" + format(
                        anc.getPosition(ax['z']) / 1000, '.15f') + "\t" + format(x, '.15f') + '\t' + format(y,
                                                                                                            '.15f') + '\t' + format(
                        r, '.15f') + '\t' + format(r_R, '.15f') + '\n')
                file.close()

    anc.close()
    #thread1.join()



def Field_scan(set_points, ramp_rate, balance_axis_index, channel_index, time_constant, balance_channel_index, filename_head):
    set_points = np.array(set_points)
    num_of_targets = len(set_points)
    print('Go to '+str(num_of_targets)+' set points')

    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)
    axis_rot = newport.NewportESP301Axis(controller,balance_axis_index-1)
    axis_rot.enable()

    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    #Connect OptiCool
    telnetObj = optc.connect_opticool()
    home_field = optc.read_field(telnetObj)[0]
    print('Current field = '+str(home_field)+' Oe')

    #Filename & title
    logbook = log_data_generator(filename_head)
    filename = logbook[0]
    savefile = logbook[1]
    if (savefile):
        file = open(filename,'a')
        file.write("Field (Oe)"+"\t"+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
        file.close()

    field_record = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

    fig = plt.figure(figsize=(8,10))
    gs = fig.add_gridspec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax1.set_xlabel('Field (Oe)')
    ax1.set_ylabel('Demod x')
    ax2.set_xlabel('Field (Oe)')
    ax2.set_ylabel('Demod y')
    ax3.set_xlabel('Field (Oe)')
    ax3.set_ylabel('R')
    draw_x, = ax1.plot([],'-o')
    draw_y, = ax2.plot([],'-o')
    draw_r, = ax3.plot([],'-o')
    fig.canvas.draw()
    fig.show()

    for i in np.arange(0,num_of_targets,1):
        if button == False:
            back_home_field = optc.set_field(telnetObj,home_field,110,0)
            break
        field = optc.read_field(telnetObj)[0]
        if (field < set_points[i]):
            ramp_up = optc.set_field(telnetObj,set_points[i],ramp_rate,0)
            time.sleep(0.5)
            #Scan
            while (field < set_points[i]-2):
                if button == False:
                    back_home_field = optc.set_field(telnetObj,home_field,110,0)
                    break
                sample0 = daq.getSample(channel_name[balance_channel_index-1] % device)
                x0 = sample0["x"][0]
                P = 1/0.3607853
                motion = -P*x0
                axis_rot.move(motion,absolute=False)
                while (axis_rot.is_motion_done==False):
                    pass
                field = optc.read_field(telnetObj)[0]
                time.sleep(time_constant*4)
                sample = daq.getSample(channel_name[channel_index-1] % device)
                sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
                x = sample["x"][0]
                y = sample["y"][0]
                r = sample["R"][0]
                field_record = np.append(field_record, field)
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_r = np.append(demod_r, r)

                #Plot
                draw_x.set_data(field_record,demod_x)
                draw_y.set_data(field_record,demod_y)
                draw_r.set_data(field_record,demod_r)
                ax1.relim()
                ax1.autoscale()
                ax2.relim()
                ax2.autoscale()
                ax3.relim()
                ax3.autoscale()
                fig.canvas.draw()
                fig.canvas.flush_events()

                #Save file
                if (savefile):
                    file = open(filename,'a')
                    file.write(format(field, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
                    file.close()

        if (field > set_points[i]):
            ramp_down = optc.set_field(telnetObj,set_points[i],ramp_rate,0)
            time.sleep(0.5)
            #Scan
            while (field > set_points[i]+2):
                if button == False:
                    back_home_field = optc.set_field(telnetObj,home_field,110,0)
                    break
                sample0 = daq.getSample(channel_name[balance_channel_index-1] % device)
                x0 = sample0["x"][0]
                P = 1/0.3607853
                motion = -P*x0
                axis_rot.move(motion,absolute=False)
                while (axis_rot.is_motion_done==False):
                    pass
                field = optc.read_field(telnetObj)[0]
                time.sleep(time_constant*4)
                sample = daq.getSample(channel_name[channel_index-1] % device)
                sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
                x = sample["x"][0]
                y = sample["y"][0]
                r = sample["R"][0]
                field_record = np.append(field_record, field)
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_r = np.append(demod_r, r)

                #Plot
                draw_x.set_data(field_record,demod_x)
                draw_y.set_data(field_record,demod_y)
                draw_r.set_data(field_record,demod_r)
                ax1.relim()
                ax1.autoscale()
                ax2.relim()
                ax2.autoscale()
                ax3.relim()
                ax3.autoscale()
                fig.canvas.draw()
                fig.canvas.flush_events()

                #Save file
                if (savefile):
                    file = open(filename,'a')
                    file.write(format(field, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
                    file.close()


    optc.disconnect_opticool(telnetObj)
    thread1.join()

def Balance_PID_single(incident_pol_angle, P, tolerance, balance_axis_index, channel_index, time_constant):
    print('Balance for', incident_pol_angle, 'incident polarization')
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    status = True
    x = 10000
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)
    axis_rot = newport.NewportESP301Axis(controller,balance_axis_index-1)
    axis_rot.enable()
    while (np.abs(x)>tolerance):
        time.sleep(time_constant*4)
        sample = daq.getSample(channel_name[channel_index-1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        print(x)
        motion = -P*x
        axis_rot.move(motion,absolute=False)
        while (axis_rot.is_motion_done==False):
            pass
    print('Balance angle = '+str(axis_rot.position))

def Balance_PID_continuous(P, balance_axis_index, channel_index, time_constant, daq, device):
    global button
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)
    axis_rot = newport.NewportESP301Axis(controller,balance_axis_index-1)
    axis_rot.enable()
    while True:
        if (button == False):
            break
        time.sleep(time_constant*4)
        sample = daq.getSample(channel_name[channel_index-1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        motion = -P*x
        axis_rot.move(motion,absolute=False)
        while (axis_rot.is_motion_done==False):
            pass

def Field_scan_PID(set_points, ramp_rate, balance_axis_index, channel_index, time_constant, P, balance_channel_index, balance_time_constant, filename_head):
    set_points = np.array(set_points)
    num_of_targets = len(set_points)
    print('Go to '+str(num_of_targets)+' set points')

    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)
    axis_rot = newport.NewportESP301Axis(controller,balance_axis_index-1)
    axis_rot.enable()

    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    #Connect OptiCool
    telnetObj = optc.connect_opticool()
    home_field = optc.read_field(telnetObj)[0]
    print('Current field = '+str(home_field)+' Oe')

    #Filename & title
    logbook = log_data_generator(filename_head)
    filename = logbook[0]
    savefile = logbook[1]
    if (savefile):
        file = open(filename,'a')
        file.write("Field (Oe)"+"\t"+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
        file.close()

    field_record = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()
    thread2 = threading.Thread(target=Balance_PID_continuous, args=(P, balance_axis_index, balance_channel_index, balance_time_constant, daq, device))
    thread2.start()

    fig = plt.figure(figsize=(8,10))
    gs = fig.add_gridspec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax1.set_xlabel('Field (Oe)')
    ax1.set_ylabel('Demod x')
    ax2.set_xlabel('Field (Oe)')
    ax2.set_ylabel('Demod y')
    ax3.set_xlabel('Field (Oe)')
    ax3.set_ylabel('R')
    draw_x, = ax1.plot([],'-o')
    draw_y, = ax2.plot([],'-o')
    draw_r, = ax3.plot([],'-o')
    fig.canvas.draw()
    fig.show()

    for i in np.arange(0,num_of_targets,1):
        field = optc.read_field(telnetObj)[0]
        if button == False:
            back_home_field = optc.set_field(telnetObj,field,110,0)
            break
        if (field < set_points[i]):
            ramp_up = optc.set_field(telnetObj,set_points[i],ramp_rate,0)
            time.sleep(0.5)
            #Scan
            while (field < set_points[i]-2):
                field = optc.read_field(telnetObj)[0]
                time.sleep(time_constant*4)
                sample = daq.getSample(channel_name[channel_index-1] % device)
                sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
                x = sample["x"][0]
                y = sample["y"][0]
                r = sample["R"][0]
                field_record = np.append(field_record, field)
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_r = np.append(demod_r, r)

                #Plot
                draw_x.set_data(field_record,demod_x)
                draw_y.set_data(field_record,demod_y)
                draw_r.set_data(field_record,demod_r)
                ax1.relim()
                ax1.autoscale()
                ax2.relim()
                ax2.autoscale()
                ax3.relim()
                ax3.autoscale()
                fig.canvas.draw()
                fig.canvas.flush_events()

                if button == False:
                    back_home_field = optc.set_field(telnetObj,field,110,0)
                    break

                #Save file
                if (savefile):
                    file = open(filename,'a')
                    file.write(format(field, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
                    file.close()

        if (field > set_points[i]):
            ramp_down = optc.set_field(telnetObj,set_points[i],ramp_rate,0)
            time.sleep(0.5)
            #Scan
            while (field > set_points[i]+2):
                field = optc.read_field(telnetObj)[0]
                time.sleep(time_constant*4)
                sample = daq.getSample(channel_name[channel_index-1] % device)
                sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
                x = sample["x"][0]
                y = sample["y"][0]
                r = sample["R"][0]
                field_record = np.append(field_record, field)
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_r = np.append(demod_r, r)

                #Plot
                draw_x.set_data(field_record,demod_x)
                draw_y.set_data(field_record,demod_y)
                draw_r.set_data(field_record,demod_r)
                ax1.relim()
                ax1.autoscale()
                ax2.relim()
                ax2.autoscale()
                ax3.relim()
                ax3.autoscale()
                fig.canvas.draw()
                fig.canvas.flush_events()

                if button == False:
                    back_home_field = optc.set_field(telnetObj,field,110,0)
                    break
                #Save file
                if (savefile):
                    file = open(filename,'a')
                    file.write(format(field, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
                    file.close()


    optc.disconnect_opticool(telnetObj)
    thread1.join()
    thread2.join()

def Vector_MOKE_Mapping(filename_head):
    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    #Filename & title
    logbook = log_data_generator(filename_head)
    filename = logbook[0]
    savefile = logbook[1]
    if (savefile):
        file = open(filename,'a')
        file.write("x coordinate (um)"+"\t"+"y coordinate (um)"+"\t"+"z coordinate (um)"+"\t"+"Mx Demod x"+"\t"+"Mx Demod y"+"\t"+"My Demod x"+"\t"+"My Demod y"+"\t"+"Mz Demod x"+"\t"+"Mz Demod y"+"\t"+"R"+"\n")
        file.close()

    #Attocube initialization
    ax = {'x':0,'y':1,'z':2}
    anc = Positioner()

    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)
    axis_rot_1 = newport.NewportESP301Axis(controller,0)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller,1)
    axis_rot_2.enable()

    #Scan parameter
    x_start = float(input('x_start_position? '))
    x_step = float(input('x_step? '))
    x_num = int(input('x_numbers? '))
    x_end = x_start + x_step * x_num

    y_start = float(input('y_start_position? '))
    y_step = float(input('y_step? '))
    y_num = int(input('y_numbers? '))
    y_end = y_start + y_step * y_num

    x_range = np.arange(x_start, x_end, x_step)
    y_range = np.arange(y_start, y_end, y_step)
    x_plot_range = np.arange(x_start, x_end + x_step, x_step)
    y_plot_range = np.arange(y_start, y_end + y_step, y_step)
    x_tor = float(input('x_tolerance? '))
    y_tor = float(input('y_tolerance? '))
    time_constant = float(input('time_constant? '))   #unit: second
    go_back = float(input('go_back? '))          #preventing hysteresis
    channel_index = int(input('channel_index? '))
    R_channel_index = int(input('R_channel_index? '))
    s_balance_angle = float(input('s_balance_angle? '))
    sp_balance_angle = float(input('sp_balance_angle? '))
    p_balance_angle = float(input('p_balance_angle? '))
    go_back_1 = float(input('rotation_go_back? '))

    scan_range_1 = np.array([0, 22.5, 45])
    scan_range_2 = np.array([s_balance_angle, sp_balance_angle, p_balance_angle])
    start_pos_1 = scan_range_1[0]
    go_back_2 = -go_back_1
    x = np.zeros(3)
    y = np.zeros(3)

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

    x_position = np.array([])
    y_position = np.array([])
    Mx_demod_x = np.array([])
    Mx_demod_y = np.array([])
    My_demod_x = np.array([])
    My_demod_y = np.array([])
    Mz_demod_x = np.array([])
    Mz_demod_y = np.array([])
    demod_x_R = np.array([])
    demod_y_R = np.array([])
    demod_r_R = np.array([])

    fig = plt.figure(figsize=(5,15))
    gs = fig.add_gridspec(4, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax4 = fig.add_subplot(gs[3, 0])
    ax1.set_xlabel('x (um)')
    ax1.set_ylabel('y (um)')
    ax1.set_xticks(x_plot_range)
    ax1.set_xticklabels(x_plot_range,rotation=90)
    ax1.set_yticks(y_plot_range)
    ax2.set_xlabel('x (um)')
    ax2.set_ylabel('y (um)')
    ax2.set_xticks(x_plot_range)
    ax2.set_xticklabels(x_plot_range,rotation=90)
    ax2.set_yticks(y_plot_range)
    ax3.set_xlabel('x (um)')
    ax3.set_ylabel('y (um)')
    ax3.set_xticks(x_plot_range)
    ax3.set_xticklabels(x_plot_range,rotation=90)
    ax3.set_yticks(y_plot_range)
    ax4.set_xlabel('x (um)')
    ax4.set_ylabel('y (um)')
    ax4.set_xticks(x_plot_range)
    ax4.set_xticklabels(x_plot_range,rotation=90)
    ax4.set_yticks(y_plot_range)

    demod_x0 = np.zeros((y_num, x_num))
    demod_y0 = np.zeros((y_num, x_num))
    demod_z0 = np.zeros((y_num, x_num))
    demod_r_R0 = np.zeros((y_num, x_num))
    x_coordinates = np.arange(x_start, x_end, x_step)
    y_coordinates = np.arange(y_start, y_end, y_step)
    X_coor, Y_coor = np.meshgrid(x_coordinates, y_coordinates)
    extent=[x_start, x_end, y_start, y_end]

    pos1 = ax1.imshow(demod_x0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos1, ax=ax1)
    pos2 = ax2.imshow(demod_y0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos2, ax=ax2)
    pos3 = ax3.imshow(demod_z0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos3, ax=ax3)
    pos4 = ax4.imshow(demod_r_R0,cmap="bwr",origin='lower',extent=extent)
    fig.colorbar(pos4, ax=ax4)
    fig.canvas.draw()
    fig.tight_layout()
    fig.show()



    for y_pos in y_range:
        if button == False:
            break
        for x_pos in x_range:
            if button == False:
                break
            if (x_pos == x_start) and (y_pos == y_start):
                x_target = x_pos-go_back
                y_target = y_pos-go_back
                anc.moveAbsolute(ax['x'], int(x_target*1000))
                anc.moveAbsolute(ax['y'], int(y_target*1000))
                x_error = np.abs(x_target-anc.getPosition(ax['x'])/1000)
                y_error = np.abs(y_target-anc.getPosition(ax['y'])/1000)
                while (x_error >= x_tor) or (y_error >= y_tor):
                    if button == False:
                        break
                    time.sleep(1)
                    #clear_output(wait=True)
                    x_error = np.abs(x_target-anc.getPosition(ax['x'])/1000)
                    y_error = np.abs(y_target-anc.getPosition(ax['y'])/1000)
                    #print(x_error)
                    #print(y_error)
                    #print(anc.getStatus(ax['x']))
                    #print(anc.getStatus(ax['y']))
                    if (x_error >= x_tor):
                        anc.moveAbsolute(ax['x'], int(x_target*1000))
                    if (y_error >= y_tor):
                        anc.moveAbsolute(ax['y'], int(y_target*1000))
                if button == False:
                    break
            anc.moveAbsolute(ax['x'], int(x_pos*1000))
            x_error = np.abs(x_pos-anc.getPosition(ax['x'])/1000)
            while (x_error >= x_tor):
                if button == False:
                    break
                time.sleep(1)
                #clear_output(wait=True)
                x_error = np.abs(x_pos-anc.getPosition(ax['x'])/1000)
                #print(x_error)
                if (x_error >= x_tor):
                    anc.moveAbsolute(ax['x'], int(x_pos*1000))
            if button == False:
                break
            anc.moveAbsolute(ax['y'], int(y_pos*1000))
            y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
            while (y_error >= y_tor):
                if button == False:
                    break
                time.sleep(1)
                #clear_output(wait=True)
                y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
                #print(y_error)
                if (y_error >= y_tor):
                    anc.moveAbsolute(ax['y'], int(y_pos*1000))
            if button == False:
                break
            print("moved to "+str(anc.getPosition(ax['x'])/1000)+","+str(anc.getPosition(ax['y'])/1000))

            for i in np.arange(0,len(scan_range_1),1):
                if button == False:
                    break
                pos_1 = scan_range_1[i]
                pos_2 = scan_range_2[i]
                if (pos_1 == start_pos_1):
                    axis_rot_1.move(pos_1-go_back_1,absolute=True)
                    axis_rot_2.move(pos_2-go_back_2,absolute=True)
                    while (axis_rot_1.is_motion_done==False) or (axis_rot_2.is_motion_done==False):
                        pass
                axis_rot_1.move(pos_1,absolute=True)
                time.sleep(0.03)
                axis_rot_2.move(pos_2,absolute=True)
                while (axis_rot_1.is_motion_done==False) or (axis_rot_2.is_motion_done==False):
                    pass
                time.sleep(time_constant*4)
                x_position = np.append(x_position, anc.getPosition(ax['x'])/1000)
                y_position = np.append(y_position, anc.getPosition(ax['y'])/1000)
                sample = daq.getSample(channel_name[channel_index-1] % device)
                sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
                x[i] = sample["x"][0]
                y[i] = sample["y"][0]

            My_x = (x[0]-x[2])/2
            Mz_x = (x[0]+x[2])/2
            Mx_x = x[1] - Mz_x
            My_y = (y[0]-y[2])/2
            Mz_y = (y[0]+y[2])/2
            Mx_y = y[1] - Mz_y
            Mx_demod_x = np.append(Mx_demod_x, Mx_x)
            My_demod_x = np.append(My_demod_x, My_x)
            Mz_demod_x = np.append(Mz_demod_x, Mz_x)
            Mx_demod_y = np.append(Mx_demod_y, Mx_y)
            My_demod_y = np.append(My_demod_y, My_y)
            Mz_demod_y = np.append(Mz_demod_y, Mz_y)

            sample_R = daq.getSample(channel_name[R_channel_index-1] % device)
            sample_R["R"] = np.abs(sample_R["x"] + 1j * sample_R["y"])
            x_R = sample_R["x"][0]
            y_R = sample_R["y"][0]
            r_R = sample_R["R"][0]
            demod_x_R = np.append(demod_x_R, x_R)
            demod_y_R = np.append(demod_y_R, y_R)
            demod_r_R = np.append(demod_r_R, r_R)

             #Plot
            length = len(Mx_demod_x)
            y_num0 = length//x_num
            x_num0 = length-y_num0*x_num
            demod_x0[:y_num0, :] = np.reshape(Mx_demod_x[:y_num0*x_num], (y_num0, x_num))
            demod_y0[:y_num0, :] = np.reshape(My_demod_x[:y_num0*x_num], (y_num0, x_num))
            demod_z0[:y_num0, :] = np.reshape(Mz_demod_x[:y_num0*x_num], (y_num0, x_num))
            demod_r_R0[:y_num0, :] = np.reshape(demod_r_R[:y_num0*x_num], (y_num0, x_num))
            if (y_num0 < y_num):
                demod_x0[y_num0, :x_num0] = Mx_demod_x[y_num0*x_num:length]
                demod_y0[y_num0, :x_num0] = My_demod_x[y_num0*x_num:length]
                demod_z0[y_num0, :x_num0] = Mz_demod_x[y_num0*x_num:length]
                demod_r_R0[y_num0, :x_num0] = demod_r_R[y_num0*x_num:length]

            pos1.set_data(demod_x0)
            pos1.set_clim(vmin = demod_x0.min(), vmax = demod_x0.max())
            pos2.set_data(demod_y0)
            pos2.set_clim(vmin = demod_y0.min(), vmax = demod_y0.max())
            pos3.set_data(demod_z0)
            pos3.set_clim(vmin = demod_z0.min(), vmax = demod_z0.max())
            pos4.set_data(demod_r_R0)
            pos4.set_clim(vmin = demod_r_R0.min(), vmax = demod_r_R0.max())
            fig.canvas.draw()
            fig.canvas.flush_events()


            if (savefile):
                file = open(filename,'a')
                file.write(format(x_position[len(x_position)-1], '.15f')+"\t"+format(y_position[len(y_position)-1], '.15f')+"\t"+format(anc.getPosition(ax['z'])/1000, '.15f')+"\t"+format(Mx_x, '.15f')+'\t'+format(Mx_y, '.15f')+'\t'+format(My_x, '.15f')+'\t'+format(My_y, '.15f')+'\t'+format(Mz_x, '.15f')+'\t'+format(Mz_y, '.15f')+'\t'+format(r_R, '.15f')+'\n')
                file.close()

    anc.close()
    thread1.join()

def initialization_lakeshore335(baud):
    return lakeshore.model_335.Model335(baud)

def read_temperature():
    lsobj = ls.initialization_lakeshore335()
    temp = ls.read_temperature(lsobj)
    ls.close_lakeshore335(lsobj)
    return float(temp)

def set_setpnttemp(Ti):
    lsobj = ls.initialization_lakeshore335()
    T = float(ls.read_temperature(lsobj))

    ls.set_ramp(lsobj, 1, 1, 0)
    time.sleep(1)
    ls.set_setpoint(lsobj, 1, Ti)
    time.sleep(1)
    T = float(ls.read_temperature(lsobj))
    #while abs(float(T) - float(Ti)) > tolerance:
    #    time.sleep(10)
    #    T = float(ls.read_temperature(lsobj))
    #    print(T)
    print(T)
    ls.close_lakeshore335(lsobj)

def set_temperature(Ti, tolerance):
    lsobj = ls.initialization_lakeshore335()
    T = float(ls.read_temperature(lsobj))

    ls.set_ramp(lsobj, 1, 1, 0)
    time.sleep(1)
    ls.set_setpoint(lsobj, 1, Ti)
    time.sleep(1)
    T = float(ls.read_temperature(lsobj))
    while abs(float(T) - float(Ti)) > tolerance:
        time.sleep(10)
        T = float(ls.read_temperature(lsobj))
        print(T)
    print(T)
    ls.close_lakeshore335(lsobj)


def TempRamp1(start_temp, end_temp, step_size, channel_index, R_channel_index, time_constant, showplot, filename_head, filename):
    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    # Attocube initialization
    ax = {'x':0, 'y':1, 'z':2}
    anc = Positioner()

    # Lakeshore initialization
    import OrensteinLab.Instrument.Lakeshore.Lakeshore335 as ls
    lsobj = ls.initialization_lakeshore335()
    ls.set_ramp(lsobj, 1, 0, 0)

    # Set time constant as specified
    daq.setDouble('/%s/demods/0/timeconstant' % device, time_constant)

    # List of temperatures to measure
    Trange = get_temp_values(start_temp, end_temp, step_size)

    # Filename & title
    totfilename=add_unique_postfix(filename_head + '\\' + filename + '.dat')
    # totfilename = filename_head + '\\' + filename + '.dat'
    file = open(totfilename, 'a')
    file = open(totfilename,'a')
    file.write("Temperature (K)"+'\t'+"x coordinate (um)"+'\t'+"y coordinate (um)"+'\t'+"z coordinate (um)"+'\t'+"Demod x"+'\t'+"Demod y"+'\t'+"R"+"Demod x (R_channel)"+'\t'+"Demod y (R_channel)"+'\t'+"R (R_channel)"+'\n')
    file.close()

    global button
    button = True
    # thread1 = threading.Thread(target=get_input)
    # thread1.start()

    temperature = np.array([])
    x_position = np.array([])
    y_position = np.array([])
    z_position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])
    demod_x_R = np.array([])
    demod_y_R = np.array([])
    demod_r_R = np.array([])

    if showplot == True:
        fig = plt.figure(figsize=(8,10))
        gs = fig.add_gridspec(3, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[2, 0])
        ax1.grid(True)
        ax2.grid(True)
        ax3.grid(True)
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('Demod x')
        ax2.set_xlabel('Temperature (K)')
        ax2.set_ylabel('Demod y')
        ax3.set_xlabel('Temperature (K)')
        ax3.set_ylabel('R')
        draw_x, = ax1.plot([],'-o')
        draw_y, = ax2.plot([],'-o')
        draw_r, = ax3.plot([],'-o')
        fig.canvas.draw()
        fig.show()

    for temp in Trange:
        if button == False:
            break
        # Change temperature
        ls.set_setpoint(lsobj, 1, temp)
        time.sleep(0.1)
        currT = []
        for m in range(60):
            if button == False:
                break
            currT.append(ls.read_temperature(lsobj))
            if m >= 2 and abs(np.mean(currT[-3:]) - temp) < 0.05:
                break
            else:
                time.sleep(1)

        temperature = np.append(temperature, ls.read_temperature(lsobj))
        x_position = np.append(x_position, anc.getPosition(ax['x'])/1000)
        y_position = np.append(y_position, anc.getPosition(ax['y'])/1000)
        z_position = np.append(z_position, anc.getPosition(ax['z'])/1000)

        sample = daq.getSample(channel_name[channel_index-1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        sample_R = daq.getSample(channel_name[R_channel_index-1] % device)
        sample_R["R"] = np.abs(sample_R["x"] + 1j * sample_R["y"])
        x_R = sample_R["x"][0]
        y_R = sample_R["y"][0]
        r_R = sample_R["R"][0]
        demod_x_R = np.append(demod_x_R, x_R)
        demod_y_R = np.append(demod_y_R, y_R)
        demod_r_R = np.append(demod_r_R, r_R)

        # Plot
        if showplot == True:
            draw_x.set_data(temperature,demod_x)
            draw_y.set_data(temperature,demod_y)
            draw_r.set_data(temperature,demod_r)
            ax1.relim()
            ax1.autoscale()
            ax2.relim()
            ax2.autoscale()
            ax3.relim()
            ax3.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

        file = open(totfilename,'a')
        file.write(format(temperature[-1], '.3f')+'\t'+format(x_position[-1], '.15f')+'\t'+format(y_position[-1], '.15f')+'\t'+format(z_position[-1], '.15f')+'\t'+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\t'+format(x_R, '.15f')+'\t'+format(y_R, '.15f')+'\t'+format(r_R, '.15f')+'\n')
        file.close()

    print("Scan finished!")
    ls.close_lakeshore335(lsobj)
    anc.close()
    button = False

def TempRamp(Ti, Tf, rate, cooldown,showplot, filename_head, filename):

    filename = log_data(filename_head, filename)
    savefile = 1


    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)



    # Lakeshore initialization
    lsobj = ls.initialization_lakeshore335()
    ls.set_ramp(lsobj, 1, 0, 0)




    file = open(filename,'a')
    #file.write('Time(s)'+'\t'+'T'+'\t'+'X'+'\t'+'Y'+'\t'+'Z'+'\t'+'Rx'+'\t'+'Ry'+'Ry'+'\n')
    file.write('Time(s)'+'\t'+'T'+'\t'+'Rx'+'\t'+'Ry'+'\t'+'R'+'\t'+'theta'+'\t'+"X"+"\t"+"Y"+"\t"+"Z"+'\t'+'I'+'\n')
    file.close()

    timeArray = np.array([])
    temperature = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])
    demod_theta = np.array([])
    demod_I = np.array([])
    posX=np.array([])
    posY=np.array([])
    posZ=np.array([])

    if showplot == True:
        fig = plt.figure(figsize=(8,10))
        gs = fig.add_gridspec(3, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[2, 0])
        ax1.grid(True)
        ax2.grid(True)
        ax3.grid(True)
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('Demod x')
        ax2.set_xlabel('Temperature (K)')
        ax2.set_ylabel('Demod y')
        ax3.set_xlabel('Temperature (K)')
        ax3.set_ylabel('R')
        draw_x, = ax1.plot([],'-o')
        draw_y, = ax2.plot([],'-o')
        draw_r, = ax3.plot([],'-o')
        fig.canvas.draw()
        fig.show()



    T=float(ls.read_temperature(lsobj))
    tolerance=0.05;

    #go to initia point

    ls.set_ramp(lsobj, 1, 1, 0)
    ls.set_setpoint(lsobj, 1, Ti)
    time.sleep(1)
    ls.set_setpoint(lsobj, 1, Ti)

    time.sleep(1)
    T=float(ls.read_temperature(lsobj))
    while abs(float(T)-float(Ti))>tolerance:
        time.sleep(10)
        T=float(ls.read_temperature(lsobj))
        print(T)

    #start ramp
    ls.set_ramp(lsobj, 1, 1, rate)
    time.sleep(1)
    ls.set_setpoint(lsobj, 1, Tf)
    time.sleep(1)

    t0=time.time()
    while  abs(float(T)-float(Tf))>tolerance:

        T=float(ls.read_temperature(lsobj))
        sample = daq.getSample('/%s/demods/0/sample' % device)
        sample2 = daq.getSample('/%s/demods/2/sample' % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        sample["theta"] = np.angle(sample["x"] + 1j * sample["y"])/np.pi*180
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        I = sample2["x"][0]
        theta = sample["theta"][0]
        temperature = np.append(temperature, T)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)
        demod_theta = np.append(demod_theta, theta)

        ts=time.time()-t0

        #posx=read_attocube('x')
        #posy=read_attocube('y')
        #posz=read_attocube('z')
        #save
        file = open(filename,'a')
        # file.write(format(ts, '.15f')+"\t"+format(T, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\t'+format(theta, '.15f')+"\t"+format(posx, '.15f')+"\t"+format(posy, '.15f')+"\t"+format(posz, '.15f')+"\t"+format(I, '.15f')+'\n')
        file.write(format(ts, '.15f')+"\t"+format(T, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\t'+format(theta, '.15f')+"\t"+format(I, '.15f')+'\n')
        file.close()


        # Plot
        if showplot == True:
            draw_x.set_data(temperature,demod_x)
            draw_y.set_data(temperature,demod_y)
            draw_r.set_data(temperature,demod_r)
            ax1.relim()
            ax1.autoscale()
            ax2.relim()
            ax2.autoscale()
            ax3.relim()
            ax3.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()
        time.sleep(1)

    if cooldown:
        ls.set_ramp(lsobj, 1, 1, 0)
        time.sleep(1)
        ls.set_setpoint(lsobj, 1, Ti)
    ls.close_lakeshore335(lsobj)

def move_attocube(axis,position, tolerance, printout):
    ax = {'x': 0, 'y': 1, 'z': 2};
    anc = Positioner();
    position=position;
    tolerance=tolerance;

    anc.moveAbsolute(ax[axis], int(position * 1000))
    error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
    while (error >= tolerance) :
        time.sleep(1)
        error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
        anc.moveAbsolute(ax[axis], int(position * 1000))
    if printout:
        print('axis ' + axis + ' arrived at', anc.getPosition(ax[axis])/1000)
    anc.close()

def read_attocube(axis):
    ax = {'x': 0, 'y': 1, 'z': 2};
    anc = Positioner();
    time.sleep(0.1)
    pos= anc.getPosition(ax[axis])/1000;
    anc.close()
    return pos

def read_allattocubes():
    ax = {'x': 0, 'y': 1, 'z': 2};
    anc = Positioner();
    time.sleep(0.1)
    pos= [anc.getPosition(0)/1000, anc.getPosition(1)/1000, anc.getPosition(2)/1000];
    anc.close()
    return pos

def read_field():
    telnetObj = optc.connect_opticool();
    real_field = optc.read_field(telnetObj)[0]
    optc.disconnect_opticool(telnetObj)
    return real_field

def reach_field(field, field_rate):
    telnetObj = optc.connect_opticool();
    ramp_up = optc.set_field(telnetObj, field, field_rate, 1)
    time.sleep(0.5)
    real_field = optc.read_field(telnetObj)[0]
    while (np.abs(field - real_field) > 1):
        time.sleep(1)
        real_field = optc.read_field(telnetObj)[0]
    optc.disconnect_opticool(telnetObj)
    return real_field

def set_field(field, field_rate):
    telnetObj = optc.connect_opticool();
    ramp_up = optc.set_field(telnetObj, field, field_rate, 1)


def MOKE_field_scan(P, channel_index, balance_channel_index, field_range, field_rate, filename_head, filename):
    filename = log_data(filename_head, filename)


    #optiCool

    telnetObj = optc.connect_opticool();
    # ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    # Attocube initialization
    ax = {'x': 0, 'y': 1, 'z': 2}
    anc = Positioner()

    axis_rot_1 = newport.NewportESP301Axis(controller, 0)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller, 1)
    axis_rot_2.enable()


    file = open(filename, 'a')
    file.write("Time [s] "+'\t'+"Field (Oe)"+'\t'+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
    file.close()
    waittime = wait_time()
    t0 = time.time()
    for field in field_range:
        # field_rate = np.abs(field-last)/10

        ramp_up = optc.set_field(telnetObj, field, field_rate, 1)
        time.sleep(0.5)
        real_field = optc.read_field(telnetObj)[0]
        while (np.abs(field - real_field) > 1):
            time.sleep(1)
            real_field = optc.read_field(telnetObj)[0]

        print('H =', real_field)

        #autobalance
        autoBalance(balance_channel_index, P);
        print('Balance angle =', float(axis_rot_2.position))

        time.sleep(0.5)

        # Quantities to be measured
        position = np.array([])
        demod_x = np.array([])
        demod_y = np.array([])
        demod_r = np.array([])

        # Scan

        time.sleep(waittime)
        sample = daq.getSample(channel_name[channel_index - 1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)
        ts = time.time()-t0
        # Save file
        file = open(filename, 'a')
        file.write(
            format(float(ts), '.15f')+ "\t"  +format(float(real_field), '.15f') + "\t" + format(x, '.15f') + '\t' + format(y, '.15f') + '\t' + format(
                r, '.15f') + '\n')
        file.close()
    optc.disconnect_opticool(telnetObj)
    anc.close()

def autoBalance(balance_channel_index, P):
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)
    x = 10000
    waittime = wait_time()

    axis_rot_1 = newport.NewportESP301Axis(controller, 0)
    axis_rot_1.enable()
    axis_rot_2 = newport.NewportESP301Axis(controller, 1)
    axis_rot_2.enable()

    while (np.abs(x) > 0.001):

        time.sleep(waittime)
        sample = daq.getSample(channel_name[balance_channel_index - 1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        motion = -P * x
        axis_rot_2.move(motion, absolute=False)
        while True:
            time.sleep(0.03)
            try:
                motion_done = axis_rot_2.is_motion_done
                break
            except ValueError:
                pass
        while (motion_done == False):
            while True:
                time.sleep(0.03)
                try:
                    motion_done = axis_rot_2.is_motion_done
                    break
                except ValueError:
                    pass


def z_scan(z_start, z_end, step_size, z_tol, channel_index, showplot, filename_head, filename):
    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)
    channel_index = 1

    # Attocube initialization
    ax = {'x': 0, 'y': 1, 'z': 2}
    anc = Positioner()
    go_back = 10

    z_range = np.arange(z_start, z_end, step_size)

    append_filename = filename + '_' + str(z_start) + '-' + str(z_end) + 'um_' + str(step_size) + 'um'
    totfilename=add_unique_postfix(filename_head + '\\' + append_filename + '.dat')
    # totfilename = filename_head + '\\' + filename + '.dat'
    file = open(totfilename, 'a')
    file.write('x: '+'\t'+str(anc.getPosition(ax['x'])/1000)+'\n')
    file.write('y: '+'\t'+str(anc.getPosition(ax['y'])/1000)+'\n')
    file.write("z position (um)" + "\t" + "Demod X" + "\t" + "Demod Y" + "\t" + "Demod R" + "\n")
    file.close()

    # Quantities to be measured
    z_position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])
    demod_aux1 = np.array([])
    demod_aux2 = np.array([])

    if showplot == True:
        fig = plt.figure(figsize=(8,10))
        gs = fig.add_gridspec(3, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[2, 0])
        ax1.grid(True)
        ax2.grid(True)
        ax3.grid(True)
        ax1.set_xlabel('z position (um)')
        ax1.set_ylabel('Demod x')
        ax2.set_xlabel('z position (um)')
        ax2.set_ylabel('Demod y')
        ax3.set_xlabel('z position (um)')
        ax3.set_ylabel('R')
        draw_x, = ax1.plot([],'-o')
        draw_y, = ax2.plot([],'-o')
        draw_r, = ax3.plot([],'-o')
        fig.canvas.draw()
        fig.show()

    # Preventing hysteresis
    z_target = z_start-go_back
    anc.moveAbsolute(ax['z'], int(z_target*1000))
    z_error = np.abs(z_target-anc.getPosition(ax['z'])/1000)
    while (z_error >= z_tol):
        time.sleep(0.1)
        z_error = np.abs(z_target-anc.getPosition(ax['z'])/1000)
        if (z_error >= z_tol):
            anc.moveAbsolute(ax['z'], int(z_target*1000))
    # Starting z scan
    for z_pos in z_range:
        anc.moveAbsolute(ax['z'], int(z_pos*1000))
        z_error = np.abs(z_pos-anc.getPosition(ax['z'])/1000)
        while (z_error >= z_tol):
            time.sleep(0.1)
            z_error = np.abs(z_pos-anc.getPosition(ax['z'])/1000)
            if (z_error >= z_tol):
                anc.moveAbsolute(ax['z'], int(z_pos*1000))

        z_position = np.append(z_position, anc.getPosition(ax['z'])/1000)
        sample = daq.getSample(channel_name[channel_index-1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        aux1 = sample["auxin0"][0]
        aux2 = sample["auxin1"][0]
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)
        demod_aux1 = np.append(demod_aux1, aux1)
        demod_aux2 = np.append(demod_aux2, aux2)

        # Save file
        file = open(totfilename, 'a')
        file.write(format(float(z_position[-1]), '.15f') + "\t" + format(x, '.15f') + '\t' + format(y,'.15f') + '\t' + format(r, '.15f') + '\n')
        file.close()

        # Plot
        if showplot==True:
            draw_x.set_data(z_position,demod_x)
            draw_y.set_data(z_position,demod_y)
            draw_r.set_data(z_position,demod_r)
            ax1.relim()
            ax1.autoscale()
            ax2.relim()
            ax2.autoscale()
            ax3.relim()
            ax3.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

    anc.close()
    print("Scan finished!")



def CorotMapReg(xmin, xmax, ymin, ymax,step,tolerance,num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, step_size_2, go_back_2, channel_index,time_constant, filename_head, filename):

    t0=time.time()
    printout=1

    timeArray = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    xArray=np.array([])
    yArray=np.array([])
    RArray=np.array([])
    waittime = wait_time()

    ax = {'x': 0, 'y': 1, 'z': 2};
    anc = Positioner();


    axis='x'
    position=xmin-10

    anc.moveAbsolute(ax[axis], int(position * 1000))
    error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
    while (error >= tolerance) :
        time.sleep(1)
        error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
        anc.moveAbsolute(ax[axis], int(position * 1000))
    if printout:
        print('axis ' + axis + ' arrived at', anc.getPosition(ax[axis])/1000)

    axis='y'
    position=ymin-10

    anc.moveAbsolute(ax[axis], int(position * 1000))
    error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
    while (error >= tolerance) :
        time.sleep(1)
        error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
        anc.moveAbsolute(ax[axis], int(position * 1000))
    if printout:
        print('axis ' + axis + ' arrived at', anc.getPosition(ax[axis])/1000)



    t0=time.time()

    i=0
    j=0
    #start the map
    Nx=int((xmax-xmin)/step+1)
    Ny=int((ymax-ymin)/step+1)
    xs=np.linspace(xmin, xmax, Nx)
    ys=np.linspace(ymin, ymax, Ny)
    for i in range(0, len(xs)):
        x=xs[i]
        for j in range(0, len(ys)):
            y=ys[j]
            clear_output(wait=True)
            #move
            axis='x'
            position=x

            anc.moveAbsolute(ax[axis], int(position * 1000))
            error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
            while (error >= tolerance) :
                time.sleep(1)
                error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
                anc.moveAbsolute(ax[axis], int(position * 1000))
            if printout:
                print('axis ' + axis + ' arrived at', anc.getPosition(ax[axis])/1000)

            axis='y'
            position=y

            anc.moveAbsolute(ax[axis], int(position * 1000))
            error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
            while (error >= tolerance) :
                time.sleep(1)
                error = np.abs(position - anc.getPosition(ax[axis]) / 1000)
                anc.moveAbsolute(ax[axis], int(position * 1000))
            if printout:
                print('axis ' + axis + ' arrived at', anc.getPosition(ax[axis])/1000)

            time.sleep(waittime)

            filename1=filename+'_'+str(i+1)+'_'+str(j+1)

            Corotate_measurementNoPos(num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, step_size_2, go_back_2, channel_index,time_constant,filename_head, filename1, 0)


            print (x, y)
    anc.close()
