#!/usr/bin/env python
# coding: utf-8
# @author: Yue Sun, UCB


# from OrensteinLab_git.configuration import CONFIG_DICT
import time
import numpy as np
import sys
import serial
import os
import json
from OrensteinLab_git.configuration import CONFIG_DICT

rotstage_port = CONFIG_DICT['Rot Stage COM Port']

# function to write a command to the instrument
def write_to_instrument(command,obj=None):
    obj, obj_passed = get_obj(obj)
    full_command = command + '\n'
    obj.write(full_command.encode()) # Send the command as bytes
    time.sleep(0.1) # small delay to ensure the command is sent
    if obj_passed==False:
        close_rotstage(obj)

# function to read a response from the instrument
def read_from_instrument(obj=None):
    obj, obj_passed = get_obj(obj)

    data = obj.readline() # read the response
    if obj_passed==False:
        close_rotstage(obj)
    return data.decode('utf-8').strip() #decode the bytes to string

# def query_instrument(command,obj=None):
#     obj, obj_passed = get_obj(obj)
#     full_command = command + '\r\n'
#     obj.write(full_command.encode()) # Send the command as bytes
#     time.sleep(0.1) # small delay to ensure the command is sent
#     data = obj.readline() # read the response

#     if obj_passed==False:
#         close_preamp(obj)
#     return data.decode().strip() #decode the bytes to string


def rotate(angle, obj=None): #rotate the stage by angle (deg)
    angle_current = load_position()
    obj, obj_passed = get_obj(obj)

    CHUNK_SIZE_DEG = 20 # break angle into blocks of 20-deg
    STEPS_PER_DEG = 795.237
    CHUNK_SIZE_STEP = int(CHUNK_SIZE_DEG * STEPS_PER_DEG)

    if (angle>360) or (angle<-360):
        ValueError('rotation angle out of range')
    else:
        remaining = int(angle * STEPS_PER_DEG)
        while abs(remaining) > 0:
            steps = CHUNK_SIZE_STEP if remaining > 0 else -CHUNK_SIZE_STEP
            if abs(remaining) < CHUNK_SIZE_STEP:
                steps = remaining

            command = f"{steps}"
            write_to_instrument(command,obj)
            response = read_from_instrument(obj)
            
            if "ACK" in response:
                angle_current = angle_current + steps/STEPS_PER_DEG
                save_position(angle_current)
                remaining -= steps
            else:
                print(f"Received unexpected response: '{response}'")
                # If we get nothing, we break to avoid an infinite loop
                if not response:
                    print("Error: Arduino is silent.")
                    break

    if obj_passed==False:
        close_rotstage(obj)
        return None
    else:
        return obj

def reset_zero():
    save_position(0.0)

def set_position(angle_target, obj=None):
    angle_current = load_position()
    d_angle = angle_target - angle_current
    rotate(d_angle,obj)
    angle_current = load_position()
    return angle_current





def save_position(angle):
    with open("motor_pos.json", "w") as f:
        json.dump({"last_angle": angle}, f)

def load_position():
    if os.path.exists("motor_pos.json"):
        with open("motor_pos.json", "r") as f:
            data = json.load(f)
            return data["last_angle"]
    return 0.0  # Default if no file exists

def initialize_rotstage():
    found_unit = False
    try:
        arduino = serial.Serial(rotstage_port, 9600, timeout=5)
        print("Waiting for Nano Every to initialize...")
        time.sleep(3) # Increased wait for the Every
        
        # Clear any junk in the buffer
        arduino.reset_input_buffer()
        return arduino
    except:
        raise ValueError('Cannot Find Arduino')
       

def close_rotstage(obj):
    obj.close()

def get_obj(obj):
    '''
    helper function for getting object both within another script or on the command line.
    '''
    obj_passed = True
    if obj==None:
        obj = initialize_rotstage()
        obj_passed=False
    return obj, obj_passed

# def set_curr(curr_set):
#     import OrensteinLab_git.devices.rot_stage as rs
#     import OrensteinLab_git.devices.zurich as lockin
#     import numpy as np
    
#     rot = rs.initialize_rotstage()
#     daq, device, props = lockin.initialize_zurich_lockin()
#     data = lockin.read_zurich_lockin(daq_objs = [daq, device, props], channels=[1,3])
    
#     curr_read = data[0]['Demod 3 r']
#     d_angle = 1
#     loop = 0
#     while (np.abs(curr_set-curr_read)/curr_set>0.01)and(loop < 100):
#         rs.rotate(d_angle,obj=rot)
#         data = lockin.read_zurich_lockin(daq_objs = [daq, device, props], channels=[1,3])
#         curr_new = data[0]['Demod 3 r']
#         print(curr_new)
#         d_angle_new = (curr_set-curr_new)*d_angle/(curr_new-curr_read)
        
#         d_angle = min(d_angle_new,10)
#         d_angle = max(d_angle,-10)
#         curr_read = curr_new
#         loop = loop + 1
    
#     rs.close_rotstage(rot)
#     return curr_read