#!/usr/bin/env python
# coding: utf-8

# In[ ]:


get_ipython().run_line_magic('matplotlib', 'notebook')
import zhinst.utils as ziutils
import instruments.newport as newport
import OrensteinLab.Instrument.OptiCool.OptiCool_Control as optc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import time
import datetime
from IPython.display import clear_output
import threading
import ipywidgets as widgets
from pyanc350.v2 import Positioner

device_id = 'DEV3537'
channel_name = ['/%s/demods/0/sample','/%s/demods/1/sample','/%s/demods/2/sample','/%s/demods/3/sample']

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
    

def Find_balance_angle(incident_pol_angle, axis_index, start_pos, step_size, num_of_steps, go_back, channel_index, time_constant, filename_head):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port="COM8", baud=921600)
    
    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)
    
    logbook = log_data_generator(filename_head)
    filename = logbook[0]
    savefile = logbook[1]
    print('Balance for', incident_pol_angle, 'incident polarization')
    if (savefile):
        file = open(filename,'a')
        file.write('Balance for '+str(incident_pol_angle)+' incident polarization'+'\n')
        file.write('\n')
        file.write("Angle (deg)"+'\t'+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
        file.close()

    #Scan parameter
    end_pos = start_pos + step_size * (num_of_steps-1)
    scan_range = np.arange(start_pos, start_pos + step_size * num_of_steps, step_size)
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
        pos_ref = np.linspace(start_pos, end_pos, 1000)
    
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
    
    
def Corotate_measurement(num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, step_size_2, go_back_2, channel_index, time_constant, filename_head):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port="COM8", baud=921600)

    #Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)


    logbook = log_data_generator(filename_head)
    filename = logbook[0]
    savefile = logbook[1]
    if (savefile):
        file = open(filename,'a')
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
        
        #Save file
        if (savefile):
            file = open(filename,'a')
            file.write(format(axis_rot_1.position, '.15f')+"\t"+format(axis_rot_2.position, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
            file.close()
    thread1.join()
    

def Pump_probe(axis_index, start_pos, step_size, num_of_steps, go_back, channel_index, time_constant, filename_head):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port="COM8", baud=921600)

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
    

def Mapping(filename_head):
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
    pos4 = ax4.imshow(demod_r_R0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
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
    
    
def Field_scan(set_points, ramp_rate, balance_axis_index, channel_index, time_constant, balance_channel_index, filename_head):
    set_points = np.array(set_points)
    num_of_targets = len(set_points)
    print('Go to '+str(num_of_targets)+' set points')
    
    controller = newport.NewportESP301.open_serial(port="COM8", baud=921600)
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
    controller = newport.NewportESP301.open_serial(port="COM8", baud=921600)
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
    controller = newport.NewportESP301.open_serial(port="COM8", baud=921600)
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
    
    controller = newport.NewportESP301.open_serial(port="COM8", baud=921600)
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
    thread2.join()

