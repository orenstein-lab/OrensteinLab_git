#!/usr/bin/env python
# coding: utf-8


from pyanc350.v2 import Positioner
import numpy as np
import time
from IPython.display import clear_output
import threading
import ipywidgets as widgets

ax = {'x':0,'y':1,'z':2}
#define a dict of axes to make things simpler



def get_input():
    global button
    while True:
        newValue = input('Stop? (Y/N) ')
        time.sleep(1)
        if (newValue == 'Y'):
            button = False
        if (button == False):
            break

def Read_Position(axis_index):
    anc = Positioner()
    clear_output(wait=True)
    print(anc.getPosition(ax[axis_index])/1000)
    anc.close()

def SingleStageMove(axis_index, real_target, tolerance, go_back):
    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()
    outs = widgets.Output()
    display(outs)
    with outs:
        anc = Positioner()
        clear_output(wait=True)
    target = real_target - go_back
    position = anc.getPosition(ax[axis_index])/1000
    error = np.abs(position-target)
    anc.moveAbsolute(ax[axis_index],target*1000)
    time.sleep(0.5)
    while error >= tolerance:
        if (button == False):
            anc.close()
            break
        with outs:
            clear_output(wait=True)
            position = anc.getPosition(ax[axis_index])/1000
            error = np.abs(position-target)
            if error >= tolerance:
                print('axis moving, currently at',anc.getPosition(ax[axis_index])/1000)
            else:
                print('axis arrived at',anc.getPosition(ax[axis_index])/1000)
        time.sleep(0.5)
    
    if (button == True):
        target = real_target
        position = anc.getPosition(ax[axis_index])/1000
        error = np.abs(position-target)
        anc.moveAbsolute(ax[axis_index],target*1000)
        time.sleep(0.5)
        while error >= tolerance:
            if (button == False):
                anc.close()
                break
            with outs:
                clear_output(wait=True)
                position = anc.getPosition(ax[axis_index])/1000
                error = np.abs(position-target)
                if error >= tolerance:
                    print('axis moving, currently at',anc.getPosition(ax[axis_index])/1000)
                else:
                    print('axis arrived at',anc.getPosition(ax[axis_index])/1000)
            time.sleep(0.5)
        anc.close()
        thread1.join()

def DoubleStageMove(axis_index_x, real_target_x, tolerance_x, go_back_x, axis_index_y, real_target_y, tolerance_y, go_back_y):
    global button
    button = True
    thread1 = threading.Thread(target=get_input)
    thread1.start()
    outs = widgets.Output()
    display(outs)
    with outs:
        anc = Positioner()
        clear_output(wait=True)
    target_x = real_target_x - go_back_x
    target_y = real_target_y - go_back_y
    position_x = anc.getPosition(ax[axis_index_x])/1000
    position_y = anc.getPosition(ax[axis_index_y])/1000
    error_x = np.abs(position_x-target_x)
    error_y = np.abs(position_y-target_y)
    anc.moveAbsolute(ax[axis_index_x],target_x*1000)
    anc.moveAbsolute(ax[axis_index_y],target_y*1000)
    time.sleep(0.5)
    while (error_x >= tolerance_x) or (error_y >= tolerance_y):
        if (button == False):
            anc.close()
            break
        with outs:
            clear_output(wait=True)
            position_x = anc.getPosition(ax[axis_index_x])/1000
            position_y = anc.getPosition(ax[axis_index_y])/1000
            error_x = np.abs(position_x-target_x)
            error_y = np.abs(position_y-target_y)
            if (error_x >= tolerance_x) or (error_y >= tolerance_y):
                print('axis moving, currently at (',anc.getPosition(ax[axis_index_x])/1000,',',anc.getPosition(ax[axis_index_y])/1000,')')
            else:
                print('axis arrived at (',anc.getPosition(ax[axis_index_x])/1000,',',anc.getPosition(ax[axis_index_y])/1000,')')
        time.sleep(0.5)
    
    if (button == True):
        target_x = real_target_x
        target_y = real_target_y
        position_x = anc.getPosition(ax[axis_index_x])/1000
        position_y = anc.getPosition(ax[axis_index_y])/1000
        error_x = np.abs(position_x-target_x)
        error_y = np.abs(position_y-target_y)
        anc.moveAbsolute(ax[axis_index_x],target_x*1000)
        anc.moveAbsolute(ax[axis_index_y],target_y*1000)
        time.sleep(0.5)
        while (error_x >= tolerance_x) or (error_y >= tolerance_y):
            if (button == False):
                anc.close()
                break
            with outs:
                clear_output(wait=True)
                position_x = anc.getPosition(ax[axis_index_x])/1000
                position_y = anc.getPosition(ax[axis_index_y])/1000
                error_x = np.abs(position_x-target_x)
                error_y = np.abs(position_y-target_y)
                if (error_x >= tolerance_x) or (error_y >= tolerance_y):
                    print('axis moving, currently at (',anc.getPosition(ax[axis_index_x])/1000,',',anc.getPosition(ax[axis_index_y])/1000,')')
                else:
                    print('axis arrived at (',anc.getPosition(ax[axis_index_x])/1000,',',anc.getPosition(ax[axis_index_y])/1000,')')
            time.sleep(0.5)
        anc.close()
        thread1.join()




