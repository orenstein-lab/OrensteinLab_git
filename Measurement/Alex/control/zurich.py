import numpy as np
import zhinst.utils as ziutils
import zhinst.core as zicore
from OrensteinLab_git.configuration import config_dict
import time
import pickle

device_id = config_dict['Zurich Lockin ID']

######################
### Core Functions ###
######################

def read_zurich_lockin(daq_objs=None, time_constant=0.3, poll_timeout=500, channel_index=1, R_channel_index=1):
    '''
    in the future, change this to output a dictionary such that function that call this can pull out any variety of infomation. alternatively, make the subscriptions flexible (ie some list of objects we want to subscribe to with a convenient default) and then the output dictionary with names that make sense. I'll then have to change upstream functions to pull out the right values, or really they should just save everything that comes out. This is a fairly straight forwward thing to implement.

    Functions I'll have to modify: autobalance, find_balance_angle, lockin_time_series, corotate_scan, rotate_scan, motor_scan
    '''

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/demods/{channel_index-1}/timeconstant' % device, time_constant)
    daq.setDouble(f'/%s/demods/{R_channel_index-1}/timeconstant' % device, time_constant)
    poll_length = time_constant
    poll_timeout = poll_timeout # ms
    poll_flags = 0
    poll_return_flat_dict = True

    # subscribe to channels and read mfli
    time.sleep(time_constant*4)
    daq.subscribe(f'/%s/demods/{channel_index-1}/sample' % device)
    daq.subscribe(f'/%s/demods/{R_channel_index-1}/sample' % device)
    mfli_dict = daq.poll(poll_length, poll_timeout, poll_flags, poll_return_flat_dict)
    data_dict = mfli_dict[f'/%s/demods/{channel_index-1}/sample' % device]
    R_dict = mfli_dict[f'/%s/demods/{R_channel_index-1}/sample' % device]
    daq.unsubscribe('*')

    # extract data from mfli_dict
    x = data_dict['x']
    y = data_dict['y']
    r = np.abs(x + 1j*y)
    x_avg = np.mean(x)
    y_avg = np.mean(y)
    r_avg = np.mean(r)
    x_R = R_dict['x']
    y_R = R_dict['y']
    r_R = np.abs(x_R + 1j*y_R)
    x_R_avg = np.mean(x_R)
    y_R_avg = np.mean(y_R)
    r_R_avg = np.mean(r_R)

    return x_avg, y_avg, r_avg, x_R_avg, y_R_avg, r_R_avg

def set_zurich_output_amplitude(amplitude, daq_objs=None, wait_time=0, channel=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/sigouts/0/amplitudes/{channel-1}'%device, amplitude)
    time.sleep(wait_time)

def read_zurich_output_amplitude(daq_objs=None, channel=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    amplitude = daq.getDouble(f'/%s/sigouts/0/amplitudes/{channel-1}'%device)
    return amplitude

def set_zurich_frequency(freq, daq_objs=None, wait_time=0, osc=2):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/oscs/{osc-1}/freq'%device, freq)
    time.sleep(wait_time)

def read_zurich_frequency(daq_objs=None, osc=2):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    freq = daq.getDouble(f'/%s/oscs/{osc-1}/freq'%device)
    return freq

def set_zurich_aux_offset(offset, daq_objs=None, wait_time=0, channel=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'%s/auxouts/{int(channel-1)}/offset'%device, offset)
    time.sleep(wait_time)

def get_zurich_aux_offset(daq_objs=None, wait_time=0, channel=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    offset = daq.getDouble(f'/%s/auxouts/{int(channel-1)}/offset'%device)
    return offset

def initialize_zurich_lockin():
    apilevel=6
    obj = ziutils.create_api_session(device_id, apilevel)
    daq, device, props = obj
    return daq, device, props

def close_zurich_lockin(obj):
    return 0




#########################
### Wrapper functions ###
#########################

def set_zurich_aux_offset_1(offset, daq_objs=None, wait_time=0):
    set_zurich_aux_offset(offset, daq_objs, wait_time, 1)

def set_zurich_aux_offset_2(offset, daq_objs=None, wait_time=0):
    set_zurich_aux_offset(offset, daq_objs, wait_time, 2)

def get_zurich_aux_offset_1(daq_objs=None, wait_time=0):
    return get_zurich_aux_offset(daq_objs, wait_time, 1)

def get_zurich_aux_offset_2(daq_objs=None, wait_time=0):
    return get_zurich_aux_offset(daq_objs, wait_time, 2)
