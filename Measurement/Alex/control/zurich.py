import numpy as np
import zhinst.utils as ziutils
import zhinst.core as zicore
from OrensteinLab_git.configuration import config_dict
import time
import pickle
import os
import dill

device_id = config_dict['Zurich Lockin ID']
ZURICH_HANDLE_FNAME = config_dict['Zurich Handle']

######################
### Core Functions ###
######################

def read_zurich_lockin(daq_objs=None, time_constant=0.3, poll_timeout=500, channel_index=1, R_channel_index=4):
    '''
    in the future, change this to output a dictionary such that function that call this can pull out any variety of infomation. alternatively, make the subscriptions flexible (ie some list of objects we want to subscribe to with a convenient default) and then the output dictionary with names that make sense. I'll then have to change upstream functions to pull out the right values, or really they should just save everything that comes out. This is a fairly straight forwward thing to implement.

    Functions I'll have to modify: autobalance, find_balance_angle, lockin_time_series, corotate_scan, rotate_scan, motor_scan
    '''

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
        daq_objs = [daq, device, props]
    else:
        daq, device, props = daq_objs

    channels = [1,2,3,4]
    for channel in channels:
        daq.setDouble(f'/{device}/demods/{channel-1}/timeconstant', time_constant)
    poll_length = time_constant
    poll_timeout = poll_timeout # ms
    poll_flags = 0
    poll_return_flat_dict = True

    # subscribe to channels and read mfli
    time.sleep(time_constant*4)
    for channel in channels:
        daq.subscribe(f'/{device}/demods/{channel-1}/sample')

    # read lockin
    mfli_dict = daq.poll(poll_length, poll_timeout, poll_flags, poll_return_flat_dict)
    daq.unsubscribe('*')
    #print(list(mfli_dict.keys()))

    data_dict = {}
    # Give data and R channel index separtely for convenience
    x = np.mean(mfli_dict[f'/{device}/demods/{channel_index-1}/sample']['x'])
    y = np.mean(mfli_dict[f'/{device}/demods/{channel_index-1}/sample']['y'])
    x_R = np.mean(mfli_dict[f'/{device}/demods/{R_channel_index-1}/sample']['x'])
    y_R = np.mean(mfli_dict[f'/{device}/demods/{R_channel_index-1}/sample']['y'])
    data_dict[f'Demod x'] = x
    data_dict[f'Demod y'] = y
    data_dict[f'Demod r'] = np.abs(x + 1j*y)
    data_dict[f'Demod phase'] = np.arctan2(y,x)
    data_dict[f'R_x'] = x_R
    data_dict[f'R_y'] = y_R
    data_dict[f'R_r'] = np.abs(x_R + 1j*y_R)
    data_dict[f'R_phase'] = np.arctan2(y_R, x_R)

    # extract data from mfli_dict
    for channel in channels:
        x = np.mean(mfli_dict[f'/{device}/demods/{channel-1}/sample']['x'])
        y = np.mean(mfli_dict[f'/{device}/demods/{channel-1}/sample']['y'])
        autophase = read_zurich_autophase(daq_objs, channel) # shouldn't change during measurement
        data_dict[f'Demod {channel} x'] = x
        data_dict[f'Demod {channel} y'] = y
        data_dict[f'Demod {channel} r'] = np.abs(x + 1j*y)
        data_dict[f'Demod {channel} phase'] = np.arctan2(y,x)
        data_dict[f'Demod {channel} autophase'] = autophase

    # add any other data from lockin to data_dict

    return data_dict

def set_zurich_autophase(angle, daq_objs=None, wait_time=0, channel=1):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/{device}/demods/{channel-1}/phaseshift', angle)
    time.sleep(wait_time)

def set_zurich_output_amplitude(amplitude, daq_objs=None, wait_time=0, channel=1):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/sigouts/0/amplitudes/{channel-1}'%device, amplitude)
    time.sleep(wait_time)

def read_zurich_autophase(daq_objs=None, channel=1):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    angle = daq.getDouble(f'/{device}/demods/{channel-1}/phaseshift')
    return angle

def read_zurich_output_amplitude(daq_objs=None, channel=1):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    amplitude = daq.getDouble(f'/%s/sigouts/0/amplitudes/{channel-1}'%device)
    return amplitude

def set_zurich_frequency(freq, daq_objs=None, wait_time=0, osc=2):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/oscs/{osc-1}/freq'%device, freq)
    time.sleep(wait_time)

def read_zurich_frequency(daq_objs=None, osc=2):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    freq = daq.getDouble(f'/%s/oscs/{osc-1}/freq'%device)
    return freq

def set_zurich_aux_offset(offset, daq_objs=None, wait_time=0, channel=1):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'%s/auxouts/{int(channel-1)}/offset'%device, offset)
    time.sleep(wait_time)

def get_zurich_aux_offset(daq_objs=None, wait_time=0, channel=1):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    offset = daq.getDouble(f'/%s/auxouts/{int(channel-1)}/offset'%device)
    return offset

def initialize_zurich_lockin():
    apilevel=6
    objs = ziutils.create_api_session(device_id, apilevel)
    daq, device, props = objs

    '''
    # future - reuse zurich connection
    try:
        if os.path.exists(ZURICH_HANDLE_FNAME):
            with open(ZURICH_HANDLE_FNAME, 'rb') as f:
                objs = dill.loads(f)
                read_zurich_frequency(objs) # test the connection
        else:
            objs = ziutils.create_api_session(device_id, apilevel)
            daq, device, props = objs
            with open(ZURICH_HANDLE_FNAME, 'wb') as f:
                dill.dumps(objs, f)
    except Exception as e:
        #print(f'initialization error: {e}')
        os.remove(ZURICH_HANDLE_FNAME)
        objs = ziutils.create_api_session(device_id, apilevel)
        daq, device, props = objs
        with open(ZURICH_HANDLE_FNAME, 'wb') as f:
            dill.dumps(objs, f)
    '''

    return daq, device, props

def close_zurich_lockin(obj):
    return 0


###############################
### Chaning Lockin Settings ###
###############################

def set_zurich_select_signal(val, channel_index=1, daq_objs=None):
    '''
    selects signal source to be associated with demodulator on channel_index
    '''

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setInt(f'/{device}/demods/{int(channel_index)-1}/adcselect', val)

def set_zurich_acfilter(val, sigin=0, daq_objs=None):
    '''
    defines input coupling for signal on sigin
    '''

    if val==True or val==False:
        pass
    else:
        raise ValueError(f'value {val} must be a 1 or 0 or a bool')

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setInt(f'/{device}/sigins/{sigin}/ac', int(val))




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
