#!/usr/bin/env python
# coding: utf-8

import instruments.newport as newport
import zhinst.utils as ziutils
import numpy as np
import time
from IPython.display import clear_output

device_id = 'DEV3537'
channel_name = ['/%s/demods/0/sample','/%s/demods/1/sample','/%s/demods/2/sample','/%s/demods/3/sample']

#ESP301 initialization
controller = newport.NewportESP301.open_serial(port="COM8", baud=921600)

def Read_Position(axis_index):
    axis_delay_1 = newport.NewportESP301Axis(controller,axis_index-1)
    axis_delay_1.enable()
    print(axis_delay_1.position)

def SingleStageMove(axis_index, pos, go_back):
    axis_delay_1 = newport.NewportESP301Axis(controller,axis_index-1)
    axis_delay_1.enable()
    print('Target position = '+str(pos))
    print('Current position =')
    dh = display(str(axis_delay_1.position), display_id=True)
    axis_delay_1.move(pos-go_back,absolute=True)
    while (axis_delay_1.is_motion_done==False):
        dh.update(str(axis_delay_1.position))
    axis_delay_1.move(pos,absolute=True)
    while (axis_delay_1.is_motion_done==False):
        dh.update(str(axis_delay_1.position))

def Corotate(axis_index_1, pos_1, go_back_1, axis_index_2, pos_2, go_back_2):
    axis_delay_1 = newport.NewportESP301Axis(controller,axis_index_1-1)
    axis_delay_1.enable()
    axis_delay_2 = newport.NewportESP301Axis(controller,axis_index_2-1)
    axis_delay_2.enable()
    print('Target position = '+str(pos_1)+' '+str(pos_2))
    print('Current position =')
    dh = display(str(axis_delay_1.position)+' '+str(axis_delay_2.position), display_id=True)
    axis_delay_1.move(pos_1-go_back_1,absolute=True)
    axis_delay_2.move(pos_2-go_back_2,absolute=True)
    while (axis_delay_1.is_motion_done==False) or (axis_delay_2.is_motion_done==False):
        dh.update(str(axis_delay_1.position)+' '+str(axis_delay_2.position))
    axis_delay_1.move(pos_1,absolute=True)
    axis_delay_2.move(pos_2,absolute=True)
    while (axis_delay_1.is_motion_done==False) or (axis_delay_2.is_motion_done==False):
        dh.update(str(axis_delay_1.position)+' '+str(axis_delay_2.position))
        
def Balance_PID_single(incident_pol_angle, P, tolerance, balance_axis_index, channel_index, time_constant):
    print('Balance for', incident_pol_angle, 'incident polarization')
    apilevel = 6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)

    status = True
    x = 10000
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




