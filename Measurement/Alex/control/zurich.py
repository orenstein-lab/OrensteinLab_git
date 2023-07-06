import numpy as np
import zhinst.utils as ziutils
import zhinst.core as zicore
from OrensteinLab_git.configuration import config_dict
import time

device_id = config_dict['Zurich Lockin ID']

######################
### Core Functions ###
######################

def read_zurich_lockin(daq_objs=None, time_constant=0.3, poll_timeout=500, channel_index=1, R_channel_index=1):

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

def set_zurich_frequency(freq, daq_objs=None, wait_time=0, osc=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/oscs/{osc-1}/freq'%device, freq)
    time.sleep(wait_time)

def read_zurich_frequency(daq_objs=None, osc=1):

    # initialize
    if daq_objs==None:
        daq, device, props = initialize_zurich_lockin()
    else:
        daq, device, props = daq_objs

    freq = daq.getDouble(f'/%s/oscs/{osc-1}/freq'%device)
    return freq

def initialize_zurich_lockin():
    apilevel=6
    (daq, device, props) = ziutils.create_api_session(device_id, apilevel)
    return daq, device, props

def close_zurich_lockin(obj):
    return 0
