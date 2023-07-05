#!/usr/bin/env python
# coding: utf-8

# In[ ]:
get_ipython().run_line_magic('matplotlib', 'notebook')
import zhinst.utils as ziutils
from zhinst.ziPython import ziListEnum
import instruments.newport as newport
import OrensteinLab.Instrument.OptiCool.OptiCool_Control as optc
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os.path
import time
import datetime
import math
from IPython.display import clear_output
import threading
import ipywidgets as widgets
from pyanc350.v2 import Positioner

# Lakeshore intialization
import OrensteinLab.Measurement.Changmin.Lakeshore_fn as lsfn

# ESP301 initialization
import OrensteinLab.Measurement.Changmin.Newport_fn as npfn

# Lock-in intialization
f_conf = open(os.path.dirname(__file__) + r'\..\..\Configuration.txt', "r")
conf_info = f_conf.read()
conf_info_split = conf_info.split('\n')
device_id = conf_info_split[0].split('\t')[1]
apilevel = 6
(daq, device, props) = ziutils.create_api_session(device_id, apilevel)
channel_name = ['/%s/demods/0/sample', '/%s/demods/1/sample', '/%s/demods/2/sample', '/%s/demods/3/sample']

ax = {'x': 0, 'y': 1, 'z': 2}

def stab_xpos(anc, x_target, x_tol):
    anc.moveAbsolute(ax['x'], int(x_target*1000))
    time.sleep(0.1)
    x_error = np.abs(x_target - anc.getPosition(ax['x'])/1000)
    while (x_error >= x_tol):
        time.sleep(0.1)
        anc.moveAbsolute(ax['x'], int(x_target*1000))
        time.sleep(0.1)
        x_error = np.abs(x_target - anc.getPosition(ax['x'])/1000)

def stab_ypos(anc, y_target, y_tol):
    anc.moveAbsolute(ax['y'], int(y_target*1000))
    time.sleep(0.1)
    y_error = np.abs(y_target - anc.getPosition(ax['y'])/1000)
    while (y_error >= y_tol):
        time.sleep(0.1)
        anc.moveAbsolute(ax['y'], int(y_target*1000))
        time.sleep(0.1)
        y_error = np.abs(y_target - anc.getPosition(ax['y'])/1000)

def stab_zpos(anc, z_target, z_tol):
    anc.moveAbsolute(ax['z'], int(z_target*1000))
    time.sleep(0.1)
    z_error = np.abs(z_target - anc.getPosition(ax['z'])/1000)
    while (z_error >= z_tol):
        time.sleep(0.1)
        anc.moveAbsolute(ax['z'], int(z_target*1000))
        time.sleep(0.1)
        z_error = np.abs(z_target - anc.getPosition(ax['z'])/1000)

def Find_balance_angle(incident_pol_angle, axis_index, start_pos, step_size, num_of_steps, go_back, channel_index, time_constant, filename_head):
    #ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

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


def corotate(num_of_steps, axis1, start_pos_1, step_size_1, go_back_1, axis2, start_pos_2, step_size_2, go_back_2, channel_index, R_channel_index, time_constant, pplot, filename_head, filename):
    #Lock-in Amplifier initialization
    # apilevel = 6
    # (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    totfilename = filename_head + '\\' + filename + '.dat'
    file = open(totfilename,'a')
    file.write("Angle_1 (deg)"+'\t'+"Angle_2 (deg)"+'\t'+"Demod x"+'\t'+"Demod y"+'\t'+"R"+'\n')
    file.close()

    scan_range_1 = np.arange(start_pos_1, start_pos_1 + step_size_1 * num_of_steps, step_size_1)
    scan_range_2 = np.arange(start_pos_2, start_pos_2 + step_size_2 * num_of_steps, step_size_2)
    global button
    button = True
    # thread1 = threading.Thread(target=get_input)
    # thread1.start()

    #Quantities to be measured
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    if pplot == True:
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
    while npfn.isMoving(axis1) == False or npfn.isMoving(axis2) == False:
        time.sleep(0.1)

    for i, pos_1 in enumerate(scan_range_1):
        if button == False:
            break
        pos_2 = scan_range_2[i]
        if i == 0:
            npfn.stab_absolute_duo(axis1, pos_1 - go_back_1, axis2, pos_2 - go_back_2)
            time.sleep(0.1)
        npfn.stab_absolute_duo(axis1, pos_1, axis2, pos_2)
        time.sleep(time_constant*6.3)
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

        position = np.append(position, pos_1)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        #Plot
        if pplot == True:
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
        file.write(format(npfn.read(axis1),'.15f')+"\t"+format(npfn.read(axis2),'.15f')+"\t"+format(x,'.15f')+'\t'+format(y,'.15f')+'\t'+format(r,'.15f')+'\t'+format(x_R,'.15f')+'\t'+format(y_R,'.15f')+'\t'+format(r_R,'.15f')+'\n')
        file.close()
    #thread1.join()
    npfn.move_absolute(axis1, start_pos_1)
    npfn.move_absolute(axis2, start_pos_2)

def Pump_probe(axis_index, start_pos, step_size, num_of_steps, go_back, time_constant, plott, filename_head, filename):

    #Lock-in Amplifier initialization
    # apilevel = 6
    # (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    #Initialize Attocube
    #ax = {'x': 0, 'y': 1, 'z': 2}
    #anc = Positioner()
    filename_tot = filename_head + '\\' + filename + '.dat'
    file = open(filename_tot, 'a')
    file.write("Position (mm)" + "\t" + "Demod x" + '\t' + "Demod y" + '\t' + "x_R" + '\t' + "y_R" + '\t' + 'aux1' + '\t' + 'aux2' + '\n')
    file.close()

    #Scan parameter
    scan_range = np.arange(start_pos, start_pos + step_size * num_of_steps, step_size)

    global button
    button = True
    # thread1 = threading.Thread(target=get_input)
    # thread1.start()

    #Quantities to be measured
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])
    if plott:
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
            npfn.stab_absolute(axis_index, pos-go_back)
            time.sleep(0.1)
        npfn.stab_absolute(axis_index, pos)
        time.sleep(time_constant*6.3)
        sample = daq.getSample('/%s/demods/0/sample' % device)
        sample_R = daq.getSample('/%s/demods/2/sample' % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        x_R = sample_R["x"][0]
        y_R = sample_R["y"][0]
        aux1 = sample["auxin0"][0]
        aux2 = sample["auxin1"][0]

        #Plot
        if plott:
            position = np.append(position, pos)
            demod_x = np.append(demod_x, x)
            demod_y = np.append(demod_y, y)
            demod_r = np.append(demod_r, r)
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
        file = open(filename_tot,'a')
        file.write(format(pos, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(x_R, '.15f')+'\t'+format(y_R, '.15f')+'\t'+format(aux1, '.15f')+'\t'+format(aux2, '.15f')+'\n')
        file.close()
    # thread1.join()

def Pump_probe_tempdep(axis_index, start_pos1, step_size1, num_of_steps1, start_pos2, step_size2, num_of_steps2, go_back, channel_index, time_constant, filename_head, num_repeats, Trange):
    # ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    # Lakeshore initialization
    import OrensteinLab.Instrument.Lakeshore.Lakeshore335 as ls
    lsobj = ls.initialization_lakeshore335()
    ls.set_ramp(lsobj, 1, 0, 0)

    # Filename & title
    logbook = log_data_generator(filename_head)
    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

    for Tset in Trange:
        # Change temperature
        ls.set_setpoint(lsobj, 1, Tset)
        time.sleep(0.1)
        currT = []
        for m in range(60):
            currT.append(ls.read_temperature(lsobj))
            if m >= 2 and abs(np.mean(currT[-3:]) - Tset) < 0.05:
                break
            else:
                time.sleep(1)

        for j in range(num_repeats):
            # Open file to save
            filename = logbook[0].replace('.dat', '_' + str(Tset) + 'K_scan' + str(j+1) + '.dat')
            savefile = logbook[1]
            if (savefile):
                file = open(filename, 'a')
                file.write("Position (mm)" + "\t" + "Demod x" + '\t' + "Demod y" + '\t' + "R" + '\n')
                file.close()

            # Scan parameter
            scan_range1 = np.arange(start_pos1, start_pos1 + step_size1 * num_of_steps1, step_size1)
            scan_range2 = np.arange(start_pos2, start_pos2 + step_size2 * num_of_steps2, step_size2)
            scan_range = np.concatenate((scan_range1, scan_range2), axis=0)
            axis_delay = newport.NewportESP301Axis(controller, axis_index - 1)
            axis_delay.enable()

            # Quantities to be measured
            position = np.array([])
            demod_x = np.array([])
            demod_y = np.array([])
            demod_r = np.array([])
            demod_aux1 = np.array([])
            demod_aux2 = np.array([])

            # Scan
            for pos in scan_range:
                if button == False:
                    break
                if (pos == start_pos1):
                    axis_delay.move(pos - go_back, absolute=True)
                    while (axis_delay.is_motion_done == False):
                        pass
                axis_delay.move(pos, absolute=True)
                while (axis_delay.is_motion_done == False):
                    pass
                time.sleep(time_constant * 4)
                sample = daq.getSample('/%s/demods/0/sample' % device)
                sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
                x = sample["x"][0]
                y = sample["y"][0]
                r = sample["R"][0]
                aux1 = sample["auxin0"][0]
                aux2 = sample["auxin1"][0]
                position = np.append(position, pos)
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_r = np.append(demod_r, r)
                demod_aux1 = np.append(demod_aux1, aux1)
                demod_aux2 = np.append(demod_aux2, aux2)

                # Save file
                if (savefile):
                    file = open(filename, 'a')
                    file.write(format(pos, '.15f') + "\t" + format(x, '.15f') + '\t' + format(y, '.15f') + '\t' + format(r, '.15f') + '\t' + format(aux1, '.15f') + '\t' + format(aux2, '.15f') + '\n')
                    file.close()
            axis_delay.move(start_pos1 - go_back, absolute=True)
    print("Scan finished!")
    anc.close()
    ls.close_lakeshore335(lsobj)
    thread1.join()

def Pump_probe_xscan(axis_index, start_pos, step_size, num_of_steps, go_back, channel_index, time_constant, filename_head, num_repeats):
    # ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    # Filename & title
    logbook = log_data_generator(filename_head)

    x_start = float(input("Start position (um)?: "))
    x_step = float(input("Step size (um)?: "))
    xnum = int(input("Number of steps?: "))
    ax = {'x':0, 'y':1, 'z':2}
    anc = Positioner()

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

    for i in range(xnum):
        # Move attocube
        x_target = x_start + i*x_step
        anc.moveAbsolute(ax['x'], int(x_target*1000))
        print(anc.getPosition(ax['x'])/1000)

        for j in range(num_repeats):
            # Open file to save
            filename = logbook[0].replace('.dat', '_' + str(x_target) + '_scan' + str(j+1) + '.dat')
            savefile = logbook[1]
            if (savefile):
                file = open(filename, 'a')
                file.write("Position (mm)" + "\t" + "Demod x" + '\t' + "Demod y" + '\t' + "R" + '\n')
                file.close()

            # Scan parameter
            scan_range = np.arange(start_pos, start_pos + step_size * num_of_steps, step_size)
            axis_delay = newport.NewportESP301Axis(controller, axis_index - 1)
            axis_delay.enable()

            # Quantities to be measured
            position = np.array([])
            demod_x = np.array([])
            demod_y = np.array([])
            demod_r = np.array([])
            demod_aux1 = np.array([])
            demod_aux2 = np.array([])

            # Scan
            for pos in scan_range:
                if button == False:
                    break
                if (pos == start_pos):
                    axis_delay.move(pos - go_back, absolute=True)
                    while (axis_delay.is_motion_done == False):
                        pass
                axis_delay.move(pos, absolute=True)
                while (axis_delay.is_motion_done == False):
                    pass
                time.sleep(time_constant * 4)
                sample = daq.getSample('/%s/demods/0/sample' % device)
                sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
                x = sample["x"][0]
                y = sample["y"][0]
                r = sample["R"][0]
                aux1 = sample["auxin0"][0]
                aux2 = sample["auxin1"][0]
                position = np.append(position, pos)
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_r = np.append(demod_r, r)
                demod_aux1 = np.append(demod_aux1, aux1)
                demod_aux2 = np.append(demod_aux2, aux2)

                # Save file
                if (savefile):
                    file = open(filename, 'a')
                    file.write(format(pos, '.15f') + "\t" + format(x, '.15f') + '\t' + format(y, '.15f') + '\t' + format(r, '.15f') + '\t' + format(aux1, '.15f') + '\t' + format(aux2, '.15f') + '\n')
                    file.close()
            axis_delay.move(start_pos -go_back, absolute=True)
    print("Scan finished!")
    anc.close()
    thread1.join()

def Pump_probe_xscan_two(axis_index, start_pos1, step_size1, num_of_steps1, start_pos2, step_size2, num_of_steps2, go_back, channel_index, time_constant, filename_head, num_repeats):
    # ESP301 initialization
    controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)

    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    # Filename & title
    logbook = log_data_generator(filename_head)

    x_start = float(input("Start position (um)?: "))
    x_step = float(input("Step size (um)?: "))
    xnum = int(input("Number of steps?: "))
    ax = {'x':0, 'y':1, 'z':2}
    anc = Positioner()

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

    for i in range(xnum):
        # Move attocube
        x_target = x_start + i*x_step
        anc.moveAbsolute(ax['x'], int(x_target*1000))
        print(anc.getPosition(ax['x'])/1000)

        for j in range(num_repeats):
            # Open file to save
            filename = logbook[0].replace('.dat', '_' + str(x_target) + '_scan' + str(j+1) + '.dat')
            savefile = logbook[1]
            if (savefile):
                file = open(filename, 'a')
                file.write("Position (mm)" + "\t" + "Demod x" + '\t' + "Demod y" + '\t' + "R" + '\n')
                file.close()

            # Scan parameter
            scan_range1 = np.arange(start_pos1, start_pos1 + step_size1 * num_of_steps1, step_size1)
            scan_range2 = np.arange(start_pos2, start_pos2 + step_size2 * num_of_steps2, step_size2)
            scan_range = np.concatenate((scan_range1, scan_range2), axis=0)
            axis_delay = newport.NewportESP301Axis(controller, axis_index - 1)
            axis_delay.enable()

            # Quantities to be measured
            position = np.array([])
            demod_x = np.array([])
            demod_y = np.array([])
            demod_r = np.array([])
            demod_aux1 = np.array([])
            demod_aux2 = np.array([])

            # Scan
            for pos in scan_range:
                if button == False:
                    break
                if (pos == start_pos1):
                    axis_delay.move(pos - go_back, absolute=True)
                    while (axis_delay.is_motion_done == False):
                        pass
                axis_delay.move(pos, absolute=True)
                while (axis_delay.is_motion_done == False):
                    pass
                time.sleep(time_constant * 4)
                sample = daq.getSample('/%s/demods/0/sample' % device)
                sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
                x = sample["x"][0]
                y = sample["y"][0]
                r = sample["R"][0]
                aux1 = sample["auxin0"][0]
                aux2 = sample["auxin1"][0]
                position = np.append(position, pos)
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_r = np.append(demod_r, r)
                demod_aux1 = np.append(demod_aux1, aux1)
                demod_aux2 = np.append(demod_aux2, aux2)

                # Save file
                if (savefile):
                    file = open(filename, 'a')
                    file.write(format(pos, '.15f') + "\t" + format(x, '.15f') + '\t' + format(y, '.15f') + '\t' + format(r, '.15f') + '\t' + format(aux1, '.15f') + '\t' + format(aux2, '.15f') + '\n')
                    file.close()
            axis_delay.move(start_pos1 -go_back, absolute=True)
    print("Scan finished!")
    anc.close()
    thread1.join()

def monitor_power(avg_time, duration, filename_head, filename):
    filename_tot = filename_head + '\\' + filename + '.dat'

    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    global button
    button = True
    # thread1 = threading.Thread(target=get_input)
    # thread1.start()

    # Quantities to be measured
    demod_t = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_x3 = np.array([])
    demod_y3 = np.array([])
    demod_aux1 = np.array([])
    demod_aux2 = np.array([])

    ti = time.time()
    for i in range(duration):
        if button == False:
            break
        # Temporary quantities within the while loop
        thisaux1 = np.array([])
        thisaux2 = np.array([])
        while True:
            dt = time.time() - ti
            sample = daq.getSample('/%s/demods/0/sample' % device)
            sample3 = daq.getSample('/%s/demods/2/sample' % device)
            if dt > (i+1)*avg_time:
                x = sample["x"][0]
                y = sample["y"][0]
                x3 = sample3["x"][0]
                y3 = sample3["y"][0]
                aux1 = thisaux1.mean()
                aux2 = thisaux2.mean()
                demod_t = np.append(demod_t, dt)
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_x3 = np.append(demod_x3, x3)
                demod_y3 = np.append(demod_y3, y3)
                demod_aux1 = np.append(demod_aux1, aux1)
                demod_aux2 = np.append(demod_aux2, aux2)
                break
            thisaux1 = np.append(thisaux1, sample["auxin0"][0])
            thisaux2 = np.append(thisaux2, sample["auxin1"][0])

        file = open(filename_tot, 'a')
        file.write(format(dt, '.15f') + "\t" + format(x, '.15f') + '\t' + format(y, '.15f') + '\t' + format(x3, '.15f') + '\t' + format(y3, '.15f') + '\t' +format(aux1, '.15f') + '\t' + format(aux2, '.15f') + '\n')
        file.close()
        # time.sleep(1 - ((time.time() - ti) % 1))

    print("Scan finished!")
    # thread1.join()

def galvo_scan(x_start, x_step, x_num, axis_index, time_constant, channel_index, channel_index3, plott, filename_head, filename):
    # Filename & title
    totfilename = filename_head + '\\' + filename + '.dat'
    file = open(totfilename, 'a')
    file.write("x (V)" + "\t" + "Demod x" + "\t" + "Demod y" + "\t" + "Demod x3" + "\t" + "Demod y3" + "\n")
    file.close()

    x_end = x_start + x_step * (x_num - 1)
    x_range = np.linspace(x_start, x_end, x_num)

    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_x3 = np.array([])
    demod_y3 = np.array([])

    if plott:
        # Prepare figure
        fig = plt.figure(figsize=(8, 3))
        gs = fig.add_gridspec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.grid(True)
        ax1.set_xlabel('X (V)')
        ax1.set_ylabel('Signal (V)')
        draw_x, = ax1.plot([], '-o')

        fig.canvas.draw()
        fig.show()
        fig.tight_layout()

    for x_num0, x_pos in enumerate(x_range):
        daq.set([["%s/AUXOUTS/%d/OFFSET" % (device_id, axis_index-1), x_pos]])
        time.sleep(0.01+ 6.3 * time_constant)
        sample = daq.getSample(channel_name[channel_index - 1] % device)
        x = sample["x"][0]
        y = sample["y"][0]
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)

        sample3 = daq.getSample(channel_name[channel_index3 - 1] % device)
        x3 = sample3["x"][0]
        y3 = sample3["y"][0]
        demod_x3 = np.append(demod_x3, x3)
        demod_y3 = np.append(demod_y3, y3)

        # Plot

        pos = daq.getDouble("/%s/AUXOUTS/%d/OFFSET" % (device_id, axis_index-1))
        if plott:
            position = np.append(position, pos)
            draw_x.set_data(position, demod_x)
            ax1.relim()
            ax1.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

        file = open(totfilename, 'a')
        file.write(format(pos, '.8f') + "\t" + format(x, '.15f') + '\t' + format(y,'.15f') + "\t" + format(x3, '.15f') + '\t' + format(y3,'.15f') + '\n')
        file.close()
    print("Scan finished!")

def galvo_scan_angle(vangle, r_start, r_step, r_num, time_constant, channel_index, channel_index3, plott, filename_head, filename):
    # Filename & title
    totfilename = filename_head + '\\' + filename + '.dat'
    file = open(totfilename, 'a')
    file.write("x (V)" + "\t" + "y (V)" + "\t" + "Demod x" + "\t" + "Demod y" + "\t" + "Demod x3" + "\t" + "Demod y3" + "\n")
    file.close()

    r_end = r_start + r_step * (r_num - 1)
    r_range = np.linspace(r_start, r_end, r_num)*1e-3

    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_x3 = np.array([])
    demod_y3 = np.array([])

    if plott:
        # Prepare figure
        fig = plt.figure(figsize=(8, 3))
        gs = fig.add_gridspec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.grid(True)
        ax1.set_xlabel('X (V)')
        ax1.set_ylabel('Signal (V)')
        draw_x, = ax1.plot([], '-o')

        fig.canvas.draw()
        fig.show()
        fig.tight_layout()

    for r_num0, r_pos in enumerate(r_range):
        daq.set([["%s/AUXOUTS/%d/OFFSET" % (device_id, 0), r_pos*np.sin(np.radians(vangle))]])
        daq.set([["%s/AUXOUTS/%d/OFFSET" % (device_id, 1), r_pos*np.cos(np.radians(vangle))]])
        time.sleep(0.01+ 6.3 * time_constant)
        sample = daq.getSample(channel_name[channel_index - 1] % device)
        x = sample["x"][0]
        y = sample["y"][0]
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)

        sample3 = daq.getSample(channel_name[channel_index3 - 1] % device)
        x3 = sample3["x"][0]
        y3 = sample3["y"][0]
        demod_x3 = np.append(demod_x3, x3)
        demod_y3 = np.append(demod_y3, y3)

        # Plot
        posx = daq.getDouble("/%s/AUXOUTS/%d/OFFSET" % (device_id, 0))
        posy = daq.getDouble("/%s/AUXOUTS/%d/OFFSET" % (device_id, 1))
        if plott:
            position = np.append(position, r_pos)
            draw_x.set_data(position, demod_x)
            ax1.relim()
            ax1.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

        file = open(totfilename, 'a')
        file.write(format(posx, '.8f') + "\t" + format(posy, '.8f') + "\t" + format(x, '.15f') + '\t' + format(y,'.15f') + "\t" + format(x3, '.15f') + '\t' + format(y3,'.15f') + '\n')
        file.close()
    print("Scan finished!")

def galvo_mapping(x_start, x_step, x_num, y_start, y_step, y_num, time_constant, channel_index, channel_index3, refresh_pix, filename_head, filename):
    # Filename & title
    totfilename = filename_head + '\\' + filename + '.dat'
    file = open(totfilename, 'a')
    file.write("x (V)" + "\t" + "y (um)" + "\t" + "Demod x" + "\t" + "Demod y" + "\t" + "Demod x3" + "\t" +  "Demod y3" + "\n")
    file.close()

    # Scan parameter
    x_end = x_start + x_step * (x_num - 1)
    y_end = y_start + y_step * (y_num - 1)

    x_range = np.linspace(x_start, x_end, x_num)
    y_range = np.linspace(y_start, y_end, y_num)

    fig = plt.figure(figsize=(5, 5))
    gs = fig.add_gridspec(1, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlabel('x (V)')
    ax1.set_ylabel('y (V)')

    demod_x0 = np.zeros((y_num, x_num))
    extent = [x_start, x_end, y_start, y_end]

    pos1 = ax1.imshow(demod_x0, cmap="bwr", origin='lower', extent=extent)
    fig.colorbar(pos1, ax=ax1)
    fig.canvas.draw()
    fig.tight_layout()
    fig.show()

    lcount = 0
    xmin = 0; xmax = 0
    for y_num0, y_pos in enumerate(y_range):
        daq.set([["%s/AUXOUTS/1/OFFSET" % (device_id), y_pos]])
        for x_num0, x_pos in enumerate(x_range):
            daq.set([["%s/AUXOUTS/0/OFFSET" % (device_id), x_pos]])
            if x_num0 == 0:
                time.sleep(0.1)
            time.sleep(0.01 + 6.3*time_constant)
            sample = daq.getSample(channel_name[channel_index-1] % device)
            sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
            x = sample["x"][0]
            y = sample["y"][0]

            sample3 = daq.getSample(channel_name[channel_index3 - 1] % device)
            x3 = sample3["x"][0]
            y3 = sample3["y"][0]

             #Plot
            demod_x0[y_num0, x_num0] = x
            if lcount == 0:
                xmin =x; xmax = x
            elif x < xmin:
                xmin = x
            elif x > xmax:
                xmax = x
            if lcount % refresh_pix == 0 or lcount == x_num*y_num - 1:
                pos1.set_data(demod_x0)
                pos1.set_clim(vmin = xmin, vmax = xmax)
                fig.canvas.draw()
                fig.canvas.flush_events()

            x_pos_val = daq.getDouble("/%s/AUXOUTS/0/OFFSET" % (device_id))
            y_pos_val = daq.getDouble("/%s/AUXOUTS/1/OFFSET" % (device_id))

            file = open(totfilename,'a')
            file.write(format(x_pos_val, '.8f')+"\t"+format(y_pos_val, '.8f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f') + "\t" + format(x3, '.15f') + '\t' + format(y3, '.15f') + '\n')
            file.close()
            lcount += 1
    print("Scan finished!")

def galvo_circle(vradius, voffx, voffy, ang_start, ang_end, ang_step, time_constant, channel_index, channel_index3, plott, filename_head, filename):
    # Filename & title
    totfilename = filename_head + '\\' + filename + '.dat'
    file = open(totfilename, 'a')
    file.write(
        "Angle (deg)" + "\t" + "Vx (mV)" + "\t" + "Vy (mV)" + "\t" + "Demod x" + "\t" + "Demod y" + "\t" + "Demod x3" + "\t" + "Demod y3" + "\n")
    file.close()

    # Scan parameter
    ang_range = np.arange(ang_start, ang_end + ang_step, ang_step)
    angles = np.array([])
    demod_x = np.array([])

    # Prepare figure
    if plott:
        fig = plt.figure(figsize=(8, 3))
        gs = fig.add_gridspec(1, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax1.grid(True)
        ax1.set_xlabel('Angle (degrees)')
        ax1.set_ylabel('Signal (V)')
        draw_x, = ax1.plot([], '-o')

    for ang_num, ang in enumerate(ang_range):
        daq.set([["%s/AUXOUTS/%d/OFFSET" % (device_id, 0), vradius*np.sin(np.radians(ang)) + voffx]])
        daq.set([["%s/AUXOUTS/%d/OFFSET" % (device_id, 1), vradius*np.cos(np.radians(ang)) + voffy]])
        time.sleep(0.1+ 6.3 * time_constant)
        sample = daq.getSample(channel_name[channel_index - 1] % device)
        x = sample["x"][0]
        y = sample["y"][0]
        demod_x = np.append(demod_x, x)

        sample3 = daq.getSample(channel_name[channel_index3 - 1] % device)
        x3 = sample3["x"][0]
        y3 = sample3["y"][0]

        # Plot
        vx = daq.getDouble("/%s/AUXOUTS/%d/OFFSET" % (device_id, 0))*1e3
        vy = daq.getDouble("/%s/AUXOUTS/%d/OFFSET" % (device_id, 1))*1e3
        if plott:
            angles = np.append(angles, ang)
            draw_x.set_data(angles, demod_x)
            ax1.relim()
            ax1.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

        file = open(totfilename, 'a')
        file.write(format(ang, '.8f') + "\t" + format(vx, '.8f') + "\t" + format(vy, '.8f') + "\t" + format(x, '.15f') + '\t' + format(y,'.15f') + "\t" + format(x3, '.15f') + '\t' + format(y3,'.15f') + '\n')
        file.close()
    print("Scan finished!")

def mapping(x_start, x_step, x_num, x_tol, y_start, y_step, y_num, y_tol, time_constant, go_back, channel_index, R_channel_index, refresh_pix, filename_head, filename):

    #Filename & title
    totfilename = filename_head + '\\' + filename + '.dat'
    file = open(totfilename,'a')
    file.write("x (um)"+"\t"+"y (um)"+"\t"+"z (um)"+"\t"+"Demod x"+"\t"+"Demod y"+"\t"+"Rx"+"\t"+"Ry"+"\n")
    file.close()

    #Attocube initialization
    ax = {'x':0,'y':1,'z':2}
    anc = Positioner()

    #Scan parameter
    x_end = x_start + x_step * x_num
    y_end = y_start + y_step * y_num

    x_range = np.arange(x_start, x_end, x_step)
    y_range = np.arange(y_start, y_end, y_step)

    fig = plt.figure(figsize=(5,10))
    gs = fig.add_gridspec(2, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax1.set_xlabel('x (um)')
    ax1.set_ylabel('y (um)')
    ax2.set_xlabel('x (um)')
    ax2.set_ylabel('y (um)')

    demod_x0 = np.zeros((y_num, x_num))
    demod_x0_R = np.zeros((y_num, x_num))
    extent=[x_start, x_end, y_start, y_end]

    pos1 = ax1.imshow(demod_x0,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos1, ax=ax1)
    pos2 = ax2.imshow(demod_x0_R,cmap="bwr",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos2, ax=ax2)
    fig.canvas.draw()
    fig.tight_layout()
    fig.show()

    lcount = 0
    for y_num0, y_pos in enumerate(y_range):
        for x_num0, x_pos in enumerate(x_range):
            if lcount == 0:
                y_target = y_pos-go_back
                anc.moveAbsolute(ax['y'], int(y_target*1000))
                y_error = np.abs(y_target - anc.getPosition(ax['y'])/1000)
                while (y_error >= y_tol):
                    time.sleep(0.1)
                    y_error = np.abs(y_target - anc.getPosition(ax['y'])/1000)
                    if (y_error >= y_tol):
                        anc.moveAbsolute(ax['y'], int(y_target*1000))
            if x_num0 == 0:
                x_target = x_pos - go_back
                anc.moveAbsolute(ax['x'], int(x_target * 1000))
                x_error = np.abs(x_target - anc.getPosition(ax['x']) / 1000)
                while (x_error >= x_tol):
                    time.sleep(0.1)
                    x_error = np.abs(x_target - anc.getPosition(ax['x'])/1000)
                    if (x_error >= x_tol):
                        anc.moveAbsolute(ax['x'], int(x_target * 1000))

            anc.moveAbsolute(ax['x'], int(x_pos*1000))
            x_error = np.abs(x_pos-anc.getPosition(ax['x'])/1000)
            while (x_error >= x_tol):
                time.sleep(0.1)
                x_error = np.abs(x_pos - anc.getPosition(ax['x'])/1000)
                if (x_error >= x_tol):
                    anc.moveAbsolute(ax['x'], int(x_pos*1000))
            anc.moveAbsolute(ax['y'], int(y_pos*1000))
            y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
            while (y_error >= y_tol):
                time.sleep(0.1)
                y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
                if (y_error >= y_tol):
                    anc.moveAbsolute(ax['y'], int(y_pos*1000))

            time.sleep(time_constant*6.3)
            x_pos_read = anc.getPosition(ax['x'])/1000
            y_pos_read = anc.getPosition(ax['y'])/1000
            z_pos_read = anc.getPosition(ax['z'])/1000

            sample = daq.getSample(channel_name[channel_index-1] % device)
            x = sample["x"][0]
            y = sample["y"][0]

            sample_R = daq.getSample(channel_name[R_channel_index-1] % device)
            x_R = sample_R["x"][0]
            y_R = sample_R["y"][0]

             #Plot
            demod_x0[y_num0, x_num0] = x
            demod_x0_R[y_num0, x_num0] = x_R
            if lcount == 0:
                xmin = x; xmax=x
                xmin3 = x_R; xmax3 = x_R
            elif x < xmin:
                xmin = x
            elif x > xmax:
                xmax = x
            if x_R < xmin3:
                xmin3 = x_R
            elif x_R > xmax3:
                xmax3 = x_R
            if lcount % refresh_pix == 0 or lcount == x_num*y_num - 1:
                pos1.set_data(demod_x0)
                pos1.set_clim(vmin = xmin, vmax = xmax)
                pos2.set_data(demod_x0_R)
                pos2.set_clim(vmin = xmin3, vmax = xmax3)
                fig.canvas.draw()
                fig.canvas.flush_events()

            file = open(totfilename,'a')
            file.write(format(x_pos_read, '.15f')+"\t"+format(y_pos_read, '.15f')+"\t"+format(z_pos_read, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(x_R, '.15f')+'\t'+format(y_R, '.15f')+'\n')
            file.close()
            lcount += 1
    anc.close()
    print("Scan finished!")

def Mapping_polscans(x_start, x_step, x_num, x_tol, y_start, y_step, y_num, y_tol, go_back, channel_index, R_channel_index, filename_head, filename):
    num_of_steps = 19; axis_index_1 = 1; start_pos_1 = 0; step_size_1 = 5; go_back_1 = 1
    axis_index_2 = 2; start_pos_2 = 0; step_size_2 = -5; go_back_2 = 1;  time_constant = 0.1; pplot = False

    # Attocube initialization
    ax = {'x': 0, 'y': 1, 'z': 2}
    anc = Positioner()

    # Scan parameter
    x_end = x_start + x_step * x_num
    y_end = y_start + y_step * y_num

    x_range = np.arange(x_start, x_end, x_step)
    y_range = np.arange(y_start, y_end, y_step)
    x_plot_range = np.arange(x_start, x_end + x_step, x_step)
    y_plot_range = np.arange(y_start, y_end + y_step, y_step)

    # global button
    button = True
    # thread1 = threading.Thread(target=get_input)
    # thread1.start()

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
                x_error = np.abs(x_target - anc.getPosition(ax['x'])/1000)
                y_error = np.abs(y_target - anc.getPosition(ax['y'])/1000)
                while (x_error >= x_tol) or (y_error >= y_tol):
                    if button == False:
                        break
                    time.sleep(0.1)
                    x_error = np.abs(x_target - anc.getPosition(ax['x'])/1000)
                    y_error = np.abs(y_target - anc.getPosition(ax['y'])/1000)
                    if (x_error >= x_tol):
                        anc.moveAbsolute(ax['x'], int(x_target*1000))
                    if (y_error >= y_tol):
                        anc.moveAbsolute(ax['y'], int(y_target*1000))
                if button == False:
                    break
            anc.moveAbsolute(ax['x'], int(x_pos*1000))
            x_error = np.abs(x_pos-anc.getPosition(ax['x'])/1000)
            while (x_error >= x_tol):
                if button == False:
                    break
                time.sleep(0.1)
                x_error = np.abs(x_pos - anc.getPosition(ax['x'])/1000)
                if (x_error >= x_tol):
                    anc.moveAbsolute(ax['x'], int(x_pos*1000))
            if button == False:
                break
            anc.moveAbsolute(ax['y'], int(y_pos*1000))
            y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
            while (y_error >= y_tol):
                if button == False:
                    break
                time.sleep(0.1)
                y_error = np.abs(y_pos-anc.getPosition(ax['y'])/1000)
                if (y_error >= y_tol):
                    anc.moveAbsolute(ax['y'], int(y_pos*1000))
            if button == False:
                break
            time.sleep(time_constant*4)
            thisfilename = filename + '_x' + str(x_pos) + '_y' + str(y_pos)
            corotate(num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1, axis_index_2, start_pos_2, step_size_2, go_back_2, channel_index, R_channel_index, time_constant, pplot, filename_head, thisfilename)
    anc.close()
    print("Scan finished!")
    # button = False
    # thread1.join()

def Mapping_AC(filename_head):
    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    # Filename & title
    logbook = log_data_generator(filename_head)
    filename = logbook[0]
    savefile = logbook[1]
    if (savefile):
        file = open(filename, 'a')
        file.write(
            "x coordinate (um)" + "\t" + "y coordinate (um)" + "\t" + "z coordinate (um)" + "\t" + "Demod x" + "\t" + "Demod y" + "\t" + "R" + "\n")
        file.close()

    # Attocube initialization
    ax = {'x': 0, 'y': 1, 'z': 2}
    anc = Positioner()

    # Scan parameter
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
    time_constant = float(input('time_constant? '))  # unit: second
    # time_avg = float(input('averaging_time')) # unit:second
    go_back = float(input('go_back? '))  # preventing hysteresis
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
    demod_aux1 = np.array([])
    demod_aux2 = np.array([])

    x_coordinates = np.arange(x_start, x_end, x_step)
    y_coordinates = np.arange(y_start, y_end, y_step)
    X_coor, Y_coor = np.meshgrid(x_coordinates, y_coordinates)
    extent = [x_start, x_end, y_start, y_end]

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
                    time.sleep(0.1)
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
                time.sleep(0.1)
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
                time.sleep(0.1)
                # clear_output(wait=True)
                y_error = np.abs(y_pos - anc.getPosition(ax['y']) / 1000)
                # print(y_error)
                if (y_error >= y_tor):
                    anc.moveAbsolute(ax['y'], int(y_pos * 1000))
            if button == False:
                break
            # print("moved to "+str(anc.getPosition(ax['x'])/1000)+","+str(anc.getPosition(ax['y'])/1000))
            time.sleep(time_constant * 4)
            x_position = np.append(x_position, anc.getPosition(ax['x']) / 1000)
            y_position = np.append(y_position, anc.getPosition(ax['y']) / 1000)
            sample = daq.getSample(channel_name[channel_index - 1] % device)
            sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
            x = sample["x"][0]
            y = sample["y"][0]
            r = sample["R"][0]
            aux1 = sample["auxin0"][0]
            aux2 = sample["auxin1"][0]
            demod_x = np.append(demod_x, x)
            demod_y = np.append(demod_y, y)
            demod_r = np.append(demod_r, r)
            demod_r = np.append(demod_r, r)
            demod_aux1 = np.append(demod_aux1, aux1)
            demod_aux2 = np.append(demod_aux2, aux2)

            if (savefile):
                file = open(filename, 'a')
                file.write(
                    format(x_position[len(x_position) - 1], '.15f') + "\t" + format(y_position[len(y_position) - 1],
                                                                                    '.15f') + "\t" + format(
                        anc.getPosition(ax['z']) / 1000, '.15f') + "\t" + format(x, '.15f') + '\t' + format(y,
                                                                                                            '.15f') + '\t' + format(
                        r, '.15f') + '\t' + format(aux1, '.15f') + '\t' + format(aux2, '.15f') + '\n')
                file.close()
    anc.close()
    print("Scan finished!")
    thread1.join()

def xscan(x_start, x_end, x_step, x_tol, time_constant, filename_head_input, filename_input):
    # # Lock-in Amplifier initialization
    # apilevel = 6
    # (daq, device, props) = ziutils.create_api_session(device_id, apilevel)
    channel_index = 1

    channel3_index = 1

    # Attocube initialization
    ax = {'x': 0, 'y': 1, 'z': 2}
    anc = Positioner()
    go_back = 10

    filename_head = filename_head_input
    filename = filename_head + '\\' + filename_input + '.dat'
    file = open(filename, 'a')
    file.write("x_position" + "\t" + "X channel" + "\t" + "Y channel" + "\t" + "X3 channel" + "\t" + "Y3 channel" + "\t" + "Aux1 channel" + "\t" + "Aux2 channel" + "\n")
    file.close()

    x_range = np.arange(x_start, x_end, x_step)

    # global button
    # button = True
    # thread1 = threading.Thread(target=get_input)
    # thread1.start()

    # Quantities to be measured
    x_position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])
    demod_x3 = np.array([])
    demod_y3 = np.array([])
    demod_r3 = np.array([])
    demod_aux1 = np.array([])
    demod_aux2 = np.array([])

    fig = plt.figure(figsize=(8, 10))
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
    draw_x, = ax1.plot([], '-o')
    draw_y, = ax2.plot([], '-o')
    draw_r, = ax3.plot([], '-o')
    fig.canvas.draw()
    fig.show()

    for x_pos in x_range:
        # if button == False:
        #     break
        if x_pos == x_start:
            x_target = x_pos - go_back
            anc.moveAbsolute(ax['x'], int(x_target * 1000))
            time.sleep(0.1)
            x_error = np.abs(x_target - anc.getPosition(ax['x']) / 1000)
            while x_error > x_tol:
                # if button == False:
                #     break
                time.sleep(0.1)
                x_error = np.abs(x_target - anc.getPosition(ax['x']) / 1000)
                if x_error >= x_tol:
                    anc.moveAbsolute(ax['x'], int(x_target * 1000))
            # if button == False:
            #     break

        anc.moveAbsolute(ax['x'], int(x_pos * 1000))
        x_temp = np.array([])
        for j in range(100):
            # if button == False:
            #     break
            time.sleep(0.1)
            x_temp = np.append(x_temp, anc.getPosition(ax['x']) / 1000)
            if j >= 2:
                x_mean = np.mean(x_temp[-3:])
                x_error = abs(x_mean - x_pos)
                if x_error < x_tol:
                    break
                elif j % 10 == 0:
                    anc.moveAbsolute(ax['x'], int(x_pos * 1000))
        # if button == False:
        #     break

        time.sleep(time_constant * 6.3)
        x_position = np.append(x_position, anc.getPosition(ax['x']) / 1000)
        sample = daq.getSample(channel_name[channel_index - 1] % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        sample3 = daq.getSample(channel_name[channel3_index - 1] % device)
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        x3 = sample3["x"][0]
        y3 = sample3["y"][0]
        aux1 = sample["auxin0"][0]
        aux2 = sample["auxin1"][0]
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)
        demod_x3 = np.append(demod_x3, x3)
        demod_y3 = np.append(demod_y3, y3)
        demod_aux1 = np.append(demod_aux1, aux1)
        demod_aux2 = np.append(demod_aux2, aux2)

        # Plot
        draw_x.set_data(x_position, demod_x)
        draw_y.set_data(x_position, demod_y)
        draw_r.set_data(x_position, demod_r)
        ax1.relim()
        ax1.autoscale()
        ax2.relim()
        ax2.autoscale()
        ax3.relim()
        ax3.autoscale()
        fig.canvas.draw()
        fig.canvas.flush_events()

        # Save file
        file = open(filename, 'a')
        file.write(format(x_position[len(x_position) - 1], '.15f') + "\t" + format(x, '.15f') + '\t' + format(y,'.15f') + '\t' + format(x3, '.15f') + '\t' + format(y3,'.15f') + '\t' + format(aux1, '.15f') + '\t' + format(aux2, '.15f') + '\n')
        file.close()
    anc.close()
    print("Scan finished!")
    # thread1.join()

def xscan_zdep(z_start, z_end, z_step, z_tol, x_start, x_end, x_step, x_tol, time_constant, scannum, filename_head_input, filename_input):
    # Lock-in Amplifier initialization
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)
    channel_index = 1
    channel3_index = 3

    # Attocube initialization
    ax = {'x': 0, 'y': 1, 'z': 2}
    anc = Positioner()
    go_back = 10

    z_range = np.arange(z_start, z_end, z_step)
    x_range = np.arange(x_start, x_end, x_step)

    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()

    for z_pos in z_range:
        if button == False:
            break
        # if z_pos == z_start:
        #     z_target = z_pos - go_back
        #     anc.moveAbsolute(ax['z'], int(z_target * 1000))
        #     time.sleep(0.1)
        #     z_error = np.abs(z_target - anc.getPosition(ax['z']) / 1000)
        #     while z_error > z_tol:
        #         if button == False:
        #             break
        #         time.sleep(0.1)
        #         z_error = np.abs(z_target - anc.getPosition(ax['z']) / 1000)
        #         if z_error >= z_tol:
        #             anc.moveAbsolute(ax['z'], int(z_target * 1000))
        #     if button == False:
        #         break

        anc.moveAbsolute(ax['z'], int(z_pos * 1000))
        z_temp = np.array([])
        for i in range(100):
            if button == False:
                break
            time.sleep(0.1)
            z_temp = np.append(z_temp, anc.getPosition(ax['z']) / 1000)
            if i>=2 and abs(np.mean(z_temp[-3:]) - z_pos) < z_tol:
                break
            if i>=2 and abs(np.mean(z_temp[-3:]) - z_pos) >= z_tol:
                anc.moveAbsolute(ax['z'], int(z_pos * 1000))

        for i in range(scannum):
            if button == False:
                break
            filename_head = filename_head_input
            filename = filename_head + '\\' + filename_input + '_z' + str(z_pos) + '_scan' + str(i+1) + '.dat'
            file = open(filename, 'a')
            file.write("x_position" + "\t" + "X channel" + "\t" + "Y channel" + "\t" + "X3 channel" + "\t" + "Y3 channel" + "\t" + "Aux1 channel" + "\t" + "Aux2 channel" + "\n")
            file.close()

            # Quantities to be measured
            x_position = np.array([])
            demod_x = np.array([])
            demod_y = np.array([])
            demod_r = np.array([])
            demod_x3 = np.array([])
            demod_y3 = np.array([])
            demod_r3 = np.array([])
            demod_aux1 = np.array([])
            demod_aux2 = np.array([])

            for x_pos in x_range:
                if button == False:
                    break
                if x_pos == x_start:
                    x_target = x_pos - go_back
                    anc.moveAbsolute(ax['x'], int(x_target * 1000))
                    time.sleep(0.1)
                    x_error = np.abs(x_target - anc.getPosition(ax['x']) / 1000)
                    while x_error > x_tol:
                        if button == False:
                            break
                        time.sleep(0.1)
                        x_error = np.abs(x_target - anc.getPosition(ax['x']) / 1000)
                        if x_error >= x_tol:
                            anc.moveAbsolute(ax['x'], int(x_target * 1000))
                    if button == False:
                        break

                anc.moveAbsolute(ax['x'], int(x_pos * 1000))
                x_temp = np.array([])
                for j in range(100):
                    if button == False:
                        break
                    time.sleep(0.1)
                    x_temp = np.append(x_temp, anc.getPosition(ax['x']) / 1000)
                    if j >= 2:
                        x_mean = np.mean(x_temp[-3:])
                        x_error = abs(x_mean - x_pos)
                        if x_error < x_tol:
                            break
                        anc.moveAbsolute(ax['x'], int(x_pos * 1000))
                if button == False:
                    break

                time.sleep(time_constant * 4)
                x_position = np.append(x_position, x_mean)
                sample = daq.getSample(channel_name[channel_index - 1] % device)
                sample3 = daq.getSample(channel_name[channel3_index - 1] % device)
                x = sample["x"][0]
                y = sample["y"][0]
                x3 = sample3["x"][0]
                y3 = sample3["y"][0]
                aux1 = sample["auxin0"][0]
                aux2 = sample["auxin1"][0]
                demod_x = np.append(demod_x, x)
                demod_y = np.append(demod_y, y)
                demod_x3 = np.append(demod_x3, x3)
                demod_y3 = np.append(demod_y3, y3)
                demod_aux1 = np.append(demod_aux1, aux1)
                demod_aux2 = np.append(demod_aux2, aux2)

                # Save file
                file = open(filename, 'a')
                file.write(format(x_position[len(x_position) - 1], '.15f') + "\t" + format(x, '.15f') + '\t' + format(y,'.15f') + '\t' + format(x3, '.15f') + '\t' + format(y3,'.15f') + '\t' + format(aux1, '.15f') + '\t' + format(aux2, '.15f') + '\n')
                file.close()
    anc.close()
    print("Scan finished!")
    thread1.join()


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

def fieldscan(Bi, Bf, Bstep, time_constant, filename_head, filename):
    # telnetObj = optc.connect_opticool()

    # estimate balance angles
    f = open(filename_head + "\\" + 'bangle_fieldscan.dat')
    data = f.read().splitlines()

    B_all = np.array([])
    bangle_all = np.array([])

    for line in data:
        B_all = np.append(B_all, float(line.split("\t")[0]))
        bangle_all = np.append(bangle_all, float(line.split("\t")[1]))
    slope, intercept = np.polyfit(B_all, bangle_all, 1)
    f.close()

    filename_tot = filename_head + '\\' + filename + '.dat'
    file = open(filename_tot,'a')
    file.write('B (Oe)' + '\t' + 'x' + '\t' + 'y' + '\t' +"Ix" + "\t" +"Iy" + '\n')
    file.close()

    B_all = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_Ix = np.array([])

    fig = plt.figure(figsize=(8,10))
    gs = fig.add_gridspec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax1.grid(True)
    ax2.grid(True)
    ax3.grid(True)
    ax1.set_xlabel('B (Oe)')
    ax1.set_ylabel('Demod X (V)')
    ax2.set_xlabel('B (Oe)')
    ax2.set_ylabel('Demod Y (V)')
    ax3.set_xlabel('B (Oe)')
    ax3.set_ylabel('Demod R (V)')
    draw_x, = ax1.plot([], '-o')
    draw_y, = ax2.plot([], '-o')
    draw_r, = ax3.plot([], '-o')
    fig.canvas.draw()
    fig.show()

    for B in np.arange(Bi, Bf + Bstep, Bstep):
        bangle = np.polyval([slope,intercept], B)
        npfn.stab_absolute(2, bangle)
        # optc.set_field(telnetObj, B, 110, 0)
        client.set_field(B, 110, client.field.approach_mode.linear)
        time.sleep(0.5)
        # readout = optc.read_field(telnetObj)
        field, status = client.get_field()
        while abs(field - B) > 1 or "Holding" not in status:
            time.sleep(0.5)
            # readout = optc.read_field(telnetObj)
            field, status = client.get_field()
        time.sleep(6.3 * time_constant)

        # B_read = optc.read_field(telnetObj)[0]
        field, _ = client.get_field()
        sample = daq.getSample('/%s/demods/0/sample' % device)
        sample2 = daq.getSample('/%s/demods/2/sample' % device)
        x = sample["x"][0]
        y = sample["y"][0]
        Ix = sample2["x"][0]
        Iy = sample2["y"][0]

        B_all = np.append(B_all, B)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_Ix = np.append(demod_Ix, Ix)

        # plot data
        draw_x.set_data(B_all, demod_x)
        draw_y.set_data(B_all, demod_y)
        draw_r.set_data(B_all, demod_Ix)
        ax1.relim()
        ax1.autoscale()
        ax2.relim()
        ax2.autoscale()
        ax3.relim()
        ax3.autoscale()
        fig.canvas.draw()
        fig.canvas.flush_events()

        file = open(filename_tot, 'a')
        file.write(format(B_read, '.15f') + "\t" + format(x, '.15f') + '\t' + format(y, '.15f')
                    + '\t' + format(Ix, '.15f') + '\t' + format(Iy, '.15f') + '\n')
        file.close()

def fieldscan_corotate(Bi, Bf, Bstep, num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1,
                       axis_index_2, start_pos_2, step_size_2, go_back_2, time_constant, pplot, filename_head, filename):
    telnetObj = optc.connect_opticool()

    # estimate balance angles
    f = open(filename_head + "\\" + 'bangle_fieldscan.dat')
    data = f.read().splitlines()
    B_all = np.array([])
    bangle_all = np.array([])
    for line in data:
        B_all = np.append(B_all, float(line.split("\t")[0]))
        bangle_all = np.append(bangle_all, float(line.split("\t")[1]))
    slope, intercept = np.polyfit(B_all, bangle_all, 1)
    f.close()

    for B in np.arange(Bi, Bf + Bstep, Bstep):
        bangle = np.polyval([slope,intercept], B)
        optc.set_field(telnetObj, B, 110, 0)
        time.sleep(0.5)
        readout = optc.read_field(telnetObj)
        while abs(readout[0] - B) > 1 or "Holding" not in readout[1]:
            time.sleep(0.5)
            readout = optc.read_field(telnetObj)
        time.sleep(6.3 * time_constant)
        filename_current = filename + '_' + str(B) + 'Oe'
        clear_output(wait=True)
        print("Current field: " + str(B) + " Oe")
        corotate(num_of_steps, axis_index_1, start_pos_1, step_size_1, go_back_1,
                 axis_index_2, start_pos_2 + bangle, step_size_2, go_back_2, 1, 3, time_constant, pplot, filename_head, filename_current)
    print("Scan finished!")

def balance_angle(angle_start, angle_end, angle_step, time_constant, plotdata):
    b_axis = 2
    go_back = 1
    demod_x = np.array([])
    angles = np.array([])
    if plotdata:
        fig = plt.figure(figsize=(5,4))
        gs = fig.add_gridspec(1, 1)
        ax = fig.add_subplot(gs[0, 0])
        ax.grid(True)
        ax.set_xlabel('Angle (degrees)')
        ax.set_ylabel('Demod X (V)')
        draw_x, = ax.plot([], '-o')
        fig.canvas.draw()
        fig.show()

    for angle in np.arange(angle_start, angle_end + angle_step, angle_step):
        if angle == angle_start:
            npfn.stab_absolute(b_axis, angle - go_back)
        time.sleep(0.1)
        npfn.stab_absolute(b_axis, angle)
        time.sleep(6.3*time_constant)
        sample = daq.getSample(channel_name[0] % device)
        angles = np.append(angles, angle)
        demod_x = np.append(demod_x, sample["x"][0])

        if plotdata:
            draw_x.set_data(angles, demod_x)
            ax.relim()
            ax.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()

    slope, intercept = np.polyfit(angles, demod_x, 1)
    bangle = -intercept / slope
    npfn.stab_absolute(b_axis, bangle - 5*angle_step - go_back)
    for i in range(6):
        time.sleep(0.1)
        npfn.stab_absolute(b_axis, bangle - (5 - i)*angle_step)

    return -intercept/slope

def load_bangle(filename_head):
    f = open(filename_head + "\\" + 'bangle_fieldscan.dat')
    data = f.read().splitlines()

    B_all = np.array([])
    bangle_all = np.array([])

    for line in data:
        B_all = np.append(B_all, float(line.split("\t")[0]))
        bangle_all = np.append(bangle_all, float(line.split("\t")[1]))
    bangle_slope, bangle_intercept = np.polyfit(B_all, bangle_all, 1)
    f.close()
    return bangle_slope, bangle_intercept

def read_bangle(client, filename_head):
    # telnetObj = optc.connect_opticool()
    # B = optc.read_field(telnetObj)[0]
    B, _ = client.get_field()

    # if len(bangle_slope) == 0 or len(bangle_intercept) ==0:

    f = open(filename_head + "\\" + 'bangle_fieldscan.dat')
    data = f.read().splitlines()

    B_all = np.array([])
    bangle_all = np.array([])

    for line in data:
        B_all = np.append(B_all, float(line.split("\t")[0]))
        bangle_all = np.append(bangle_all, float(line.split("\t")[1]))
    bangle_slope, bangle_intercept = np.polyfit(B_all, bangle_all, 1)
    f.close()

    bangle = np.polyval([bangle_slope, bangle_intercept], B)
    npfn.stab_absolute(2, bangle)
    return bangle

def bangle_fieldscan(Bi, Bf, Bstep, angle_start, angle_end, angle_step, mult, time_constant, plotdata, client, filename_head):
    # telnetObj = optc.connect_opticool()
    for B in np.arange(Bi, Bf + Bstep, Bstep):
        # optc.set_field(telnetObj, B, 110, 0)
        client.set_field(B, 110, client.field.approach_mode.linear)
        time.sleep(0.5)
        field, _ = client.get_field()
        while abs(field - B) > 1:
            time.sleep(0.5)
            field, _ = client.get_field()
        clear_output(wait=True)
        bangle = balance_angle(angle_start + mult*B/1e4, angle_end + mult*B/1e4, angle_step, time_constant, plotdata)
        file = open(filename_head + '\\' + 'bangle_fieldscan.dat', 'a')
        file.write(format(B, '.4f') + "\t" + format(bangle, '.4f') + '\n')
        file.close()
        print(bangle)

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

def TempRamp(Ti, Tf, rate, dt, filename_head, filename, cooldown, plotdata):
    filename_tot = filename_head + '\\' + filename + '.dat'
    # file = open(filename_tot,'a')
    # file.write('Time(s)'+'\t'+'T'+'\t'+'x'+'\t'+'y'+'\t'+"Rx"+"\t"+"Ry"+'\n')
    # file.close()

    temperature = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    #go to initial point
    lsfn.set_ramp(0)
    time.sleep(0.1)
    lsfn.stab_temperature(Ti, 0.05, 3)
    time.sleep(0.1)

    #start ramp
    lsfn.set_ramp(rate)
    time.sleep(0.1)
    lsfn.set_setpoint(Tf)
    time.sleep(0.1)

    t0=time.time()
    count = 1
    T = float(lsfn.read_temperature())
    time.sleep(0.1)

    if plotdata:
        fig = plt.figure(figsize=(8,10))
        gs = fig.add_gridspec(3, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[1, 0])
        ax3 = fig.add_subplot(gs[2, 0])
        ax1.grid(True)
        ax2.grid(True)
        ax3.grid(True)
        ax1.set_xlabel('Temperature (K)')
        ax1.set_ylabel('Demod X (V)')
        ax2.set_xlabel('Temperature (K)')
        ax2.set_ylabel('Demod Y (V)')
        ax3.set_xlabel('Temperature (K)')
        ax3.set_ylabel('Demod R (V)')

        draw_x, = ax1.plot([], '-o')
        draw_y, = ax2.plot([], '-o')
        draw_r, = ax3.plot([], '-o')

        fig.canvas.draw()
        fig.show()

    while  abs(float(T)-float(Tf)) > 0.05:
        ts = time.time() - t0
        if ts > count*dt:
            sample = daq.getSample('/%s/demods/0/sample' % device)
            sample2 = daq.getSample('/%s/demods/2/sample' % device)
            sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
            # sample["theta"] = np.angle(sample["x"] + 1j * sample["y"]) / np.pi * 180
            x = sample["x"][0]
            y = sample["y"][0]
            r = sample["R"][0]
            Ix = sample2["x"][0]
            Iy = sample2["y"][0]
            # if count % 10 == 1:
            T = float(lsfn.read_temperature())
            if T >= 18 and T <=22:
                lsfn.set_ramp(rate)

            temperature = np.append(temperature, T)
            demod_x = np.append(demod_x, x)
            demod_y = np.append(demod_y, y)
            demod_r = np.append(demod_r, r)

            file = open(filename_tot, 'a')
            file.write(format(ts, '.15f') + "\t" + format(T, '.15f') + "\t" + format(x, '.15f') + '\t'
            + format(y, '.15f') + '\t'  + format(Ix, '.15f') + '\t' + format(Iy, '.15f') + '\n')
            file.close()
            count += 1

        if plotdata:
            draw_x.set_data(temperature, demod_x)
            draw_y.set_data(temperature, demod_y)
            draw_r.set_data(temperature, demod_r)
            ax1.relim()
            ax1.autoscale()
            ax2.relim()
            ax2.autoscale()
            ax3.relim()
            ax3.autoscale()
            fig.canvas.draw()
            fig.canvas.flush_events()
            #plt.show()
            #clear_output(wait=True)
            #time.sleep(2)
        #elif printt:
            #clear_output(wait=True)
            #print(T)
            #time.sleep(0.3)
        #else:
            #time.sleep(0.3)
    print('Scan finished!')
    if cooldown:
        lsfn.set_ramp(0)
        time.sleep(0.1)
        lsfn.set_setpoint(Ti)

def tempramp_cont(Ti, Tf, Trate, total_duration, module_sampling_rate, burst_duration, filename_head, filename, cooldown):
    filename_tot = filename_head + '\\' + filename + '.dat'
    file = open(filename_tot, 'a')
    file.write('Temperature (K)'+'\t'+'time (s)'+'\t'+'x1'+'\t'+'y1'+'\t'+'x3'+'\t'+'y3'+'\n')
    file.close()

    #go to initial point
    lsfn.set_ramp(0)
    time.sleep(0.1)
    lsfn.stab_temperature(Ti, 0.05, 3)
    time.sleep(0.1)

    #start temperature ramp
    lsfn.set_ramp(Trate)
    time.sleep(0.1)
    lsfn.set_setpoint(Tf)
    time.sleep(0.1)

    # The list of signal paths that we would like to record in the module.
    demod_path = f"/{device}/demods/0/sample"
    demod_path2 = f"/{device}/demods/2/sample"
    signal_paths = []
    signal_paths.append(demod_path + ".x")  # The demodulator X output.
    signal_paths.append(demod_path + ".y")  # The demodulator Y output.
    signal_paths.append(demod_path + ".auxin0")  # The demodulator aux1 output.
    # signal_paths.append(demod_path2 + ".x")  # The demodulator X output.
    # signal_paths.append(demod_path2 + ".y")  # The demodulator Y output.

    # Check the device has demodulators.
    flags = ziListEnum.recursive | ziListEnum.absolute | ziListEnum.streamingonly
    streaming_nodes = daq.listNodes(f"/{device}", flags)
    if demod_path.upper() not in streaming_nodes:
        print(
            f"Device {device} does not have demodulators. Please modify the example to specify",
            "a valid signal_path based on one or more of the following streaming nodes: ",
            "\n".join(streaming_nodes),
        )
        raise Exception("Demodulator streaming nodes unavailable - see the message above for more information.")

    # Defined the total time we would like to record data for and its sampling rate.
    # total_duration: Time in seconds: This examples stores all the acquired data in the `data` dict - remove this
    # continuous storing in read_data_update_plot before increasing the size of total_duration!
    num_cols = int(np.ceil(module_sampling_rate * burst_duration))
    num_bursts = int(np.ceil(total_duration / burst_duration))

    # Create an instance of the Data Acquisition Module.
    daq_module = daq.dataAcquisitionModule()

    # Configure the Data Acquisition Module.
    # Set the device that will be used for the trigger - this parameter must be set.
    daq_module.set("device", device)

    # Specify continuous acquisition (type=0).
    daq_module.set("type", 0)

    # 2 = Linear interpolation.
    daq_module.set("grid/mode", 2)
    daq_module.set("count", num_bursts)
    daq_module.set("duration", burst_duration)
    daq_module.set("grid/cols", num_cols)
    # daq_module.set("save/fileformat", 1)
    # daq_module.set("save/filename", filename_tot)
    # daq_module.set("save/saveonread", 1)

    data = {}
    # A dictionary to store all the acquired data.
    for signal_path in signal_paths:
        daq_module.subscribe(signal_path)
        data[signal_path] = []

    clockbase = float(daq.getInt(f"/{device}/clockbase"))
    ts0 = np.nan
    read_count = 0

    def read_data_update_plot(data, timestamp0):
        """
        Read the acquired data out from the module and plot it. Raise an
        AssertionError if no data is returned.
        """
        data_read = daq_module.read(True)
        returned_signal_paths = [signal_path.lower() for signal_path in data_read.keys()]
        progress = daq_module.progress()[0]
        tott = [[] for _ in range(len(signal_paths))]
        totvalue = [[] for _ in range(len(signal_paths))]
        # Loop over all the subscribed signals:
        for sindex, signal_path in enumerate(signal_paths):
            if signal_path.lower() in returned_signal_paths:
                # Loop over all the bursts for the subscribed signal. More than
                # one burst may be returned at a time, in particular if we call
                # read() less frequently than the burst_duration.
                for index, signal_burst in enumerate(data_read[signal_path.lower()]):
                    if np.any(np.isnan(timestamp0)):
                        # Set our first timestamp to the first timestamp we obtain.
                        timestamp0 = signal_burst["timestamp"][0, 0]
                        # print('nan!')
                    # Convert from device ticks to time in seconds.
                    t = (signal_burst["timestamp"][0, :] - timestamp0) / clockbase
                    value = signal_burst["value"][0, :]
                    # if do_plot:
                    #     ax.plot(t, value)
                    num_samples = len(signal_burst["value"][0, :])
                    dt = (signal_burst["timestamp"][0, -1] - signal_burst["timestamp"][0, 0]) / clockbase
                    data[signal_path].append(signal_burst)

                    tott[sindex] = np.append(tott[sindex], t)
                    totvalue[sindex] = np.append(totvalue[sindex],value)
                    # print(
                    #     f"Read: {read_count}, progress: {100 * progress:.2f}%.",
                    #     f"Burst {index}: {signal_path} contains {num_samples} spanning {dt:.2f} s.",
                    # )
            else:
                # Note: If we read before the next burst has finished, there may be no new data.
                # No action required.
                pass

        # Update the plot.
        # if do_plot:
        #     ax.set_title(f"Progress of data acquisition: {100 * progress:.2f}%.")
        #     plt.pause(0.01)
        #     fig.canvas.draw()
        # print(timestamp0)
        return data, timestamp0, tott, totvalue

    # Start recording data.
    daq_module.execute()

    # Record data in a loop with timeout.
    timeout = 1.5 * total_duration
    t0_measurement = time.time()
    # The maximum time to wait before reading out new data.
    t_update = 0.9 * burst_duration
    while not daq_module.finished():
        t0_loop = time.time()
        if time.time() - t0_measurement > timeout:
            raise Exception(
                f"Timeout after {timeout} s - recording not complete. Are the streaming nodes enabled? "
                "Has a valid signal_path been specified?"
            )
        T1 = lsfn.read_temperature()
        data, ts0, totdt, totvalue = read_data_update_plot(data, ts0)
        T2 = lsfn.read_temperature()
        T = (T1 + T2)/2
        # print(T)
        read_count += 1
        # save file
        f = open(filename_tot, 'a')
        np.savetxt(f, np.transpose(np.vstack((T * np.ones(len(totdt[0])), totdt[0], totvalue))), delimiter='\t', fmt='%.15f')
        f.close()
        # We don't need to update too quickly.
        time.sleep(max(0, t_update - (time.time() - t0_loop)))

    # There may be new data between the last read() and calling finished().
    # T1 = lsfn.read_temperature()
    data, _, totdt, totvalue = read_data_update_plot(data, ts0)
    # T2 = lsfn.read_temperature()
    # T = (T1 + T2)/2
    f = open(filename_tot, 'a')
    np.savetxt(f, np.transpose(np.vstack((T * np.ones(len(totdt[0])), totdt[0], totvalue))), delimiter='\t', fmt='%.15f')
    # f.write("\n")
    f.close()
    print('Scan finished! read_count = ' + str(read_count))

    if cooldown:
        lsfn.set_ramp(0)
        time.sleep(0.1)
        lsfn.set_setpoint(Ti)

    return data
