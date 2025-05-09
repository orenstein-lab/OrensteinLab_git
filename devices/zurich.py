
import zhinst.utils as ziutils
import zhinst.core as zicore
from OrensteinLab_git.configuration import CONFIG_DICT
import numpy as np
import time
import pickle
import os
import dill

DEVICE_ID = CONFIG_DICT['Zurich Lockin ID']
DEVICE_ID2 = CONFIG_DICT['Zurich Lockin ID 2']
ZURICH_HANDLE_FNAME = CONFIG_DICT['Zurich Handle']

######################
### Core Functions ### Defaults to work with lockin at DEVICE_ID
######################

def read_zurich_lockin(daq_objs=None, time_constant=0.3, poll_length=None, poll_timeout=500, wait_factor=3, channels=[1,2,3,4], device_id=DEVICE_ID, name=""):
    '''
    in the future, change this to output a dictionary such that function that call this can pull out any variety of infomation. alternatively, make the subscriptions flexible (ie some list of objects we want to subscribe to with a convenient default) and then the output dictionary with names that make sense. I'll then have to change upstream functions to pull out the right values, or really they should just save everything that comes out. This is a fairly straight forwward thing to implement.

    Functions I'll have to modify: autobalance, find_balance_angle, lockin_time_series, corotate_scan, rotate_scan, motor_scan
    '''

    # initialize
    obj_passed = True
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
        daq_objs = [daq, device, props]
        obj_passed = False
    else:
        daq, device, props = daq_objs

    for ii, channel in enumerate(channels):
        if type(time_constant) is float or type(time_constant) is int:
            daq.setDouble(f'/{device}/demods/{channel-1}/timeconstant', time_constant)
        else:
            daq.setDouble(f'/{device}/demods/{channel-1}/timeconstant', time_constant[ii])
    if poll_length==None:
        poll_length = np.max(time_constant)
    poll_timeout = poll_timeout # ms
    poll_flags = 0
    poll_return_flat_dict = True

    # subscribe to channels and read mfli
    time.sleep(np.max(time_constant)*wait_factor)
    for channel in channels:
        daq.subscribe(f'/{device}/demods/{channel-1}/sample')

    # read lockin
    daq.sync()
    mfli_dict = daq.poll(poll_length, poll_timeout, poll_flags, poll_return_flat_dict)
    daq.unsubscribe('*')
    #print(list(mfli_dict.keys()))

    # extract data from mfli_dict
    data_dict = {}
    for channel in channels:
        x = np.mean(mfli_dict[f'/{device}/demods/{channel-1}/sample']['x'])
        y = np.mean(mfli_dict[f'/{device}/demods/{channel-1}/sample']['y'])
        autophase = daq.getDouble(f'/{device}/demods/{channel-1}/phaseshift')
        osc = daq.getInt(f'/{device}/demods/{channel-1}/oscselect')
        freq = daq.getDouble(f'/{device}/oscs/{osc}/freq')
        data_dict[f'{name}Demod {channel} x'] = x
        data_dict[f'{name}Demod {channel} y'] = y
        data_dict[f'{name}Demod {channel} r'] = np.abs(x + 1j*y)
        data_dict[f'{name}Demod {channel} phase'] = np.arctan2(y,x)
        data_dict[f'{name}Demod {channel} autophase'] = autophase
        data_dict[f'{name}Demod {channel} oscillator'] = osc+1
        data_dict[f'{name}Demod {channel} frequency'] = freq

    # add any other data from lockin to data_dict
    if obj_passed==False:
        close_zurich_lockin(daq_objs)
        return data_dict, None
    else:
        return data_dict, daq_objs

def set_zurich_osc(osc, daq_objs=None, wait_time=0, channel=1, check_stability=True, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    daq.setInt(f'/{device}/demods/{channel-1}/oscselect', osc-1)
    time.sleep(wait_time)

    return daq_objs

def set_zurich_autophase(angle, daq_objs=None, wait_time=0, channel=1, check_stability=True, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/{device}/demods/{channel-1}/phaseshift', angle)
    time.sleep(wait_time)

    return daq_objs

def set_zurich_output_amplitude(amplitude, daq_objs=None, wait_time=0, channel=1, check_stability=True, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/sigouts/0/amplitudes/{channel-1}'%device, amplitude)
    time.sleep(wait_time)

    return daq_objs

def read_zurich_osc(daq_objs=None, channel=1, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    osc = daq.getInt(f'/{device}/demods/{channel-1}/oscselect')
    return osc+1, daq_objs

def read_zurich_autophase(daq_objs=None, channel=1, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    angle = daq.getDouble(f'/{device}/demods/{channel-1}/phaseshift')
    return angle, daq_objs

def read_zurich_output_amplitude(daq_objs=None, channel=1, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    amplitude = daq.getDouble(f'/%s/sigouts/0/amplitudes/{channel-1}'%device)
    return amplitude, daq_objs

def set_zurich_frequency(freq, daq_objs=None, wait_time=0, osc=2, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'/%s/oscs/{osc-1}/freq'%device, freq)
    time.sleep(wait_time)

    return daq_objs

def read_zurich_frequency(daq_objs=None, osc=2, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    freq = daq.getDouble(f'/%s/oscs/{osc-1}/freq'%device)
    return freq, daq_objs

def set_zurich_aux_offset(offset, daq_objs=None, wait_time=0, channel=1, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    daq.setDouble(f'%s/auxouts/{int(channel-1)}/offset'%device, offset)
    time.sleep(wait_time)

    return daq_objs

def get_zurich_aux_offset(daq_objs=None, wait_time=0, channel=1, device_id=DEVICE_ID):

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin(device_id)
    else:
        daq, device, props = daq_objs

    offset = daq.getDouble(f'/%s/auxouts/{int(channel-1)}/offset'%device)
    return offset, daq_objs

def set_zurich_output_ramp(v, daq_objs=None, vstep=0.5, channel=1, wait_time=2, check_stability=True, device_id=DEVICE_ID):
    '''
    For slowly stepping aux out voltage on lockin to a new setpoint.
    '''

    if daq_objs is None:
        daq_objs = initialize_zurich_lockin(device_id)
    vcurr, daq_objs = read_zurich_output_amplitude(daq_objs=daq_objs, channel=channel, device_id=device_id)
    vcurr = round(vcurr,1)
    sgn = np.sign(v - vcurr)
    if vcurr==v:
        pass
    else:
        voltages = np.arange(vcurr,v+sgn*vstep, sgn*vstep)
        for vv in voltages:
            set_zurich_output_amplitude(vv,wait_time=wait_time, daq_objs=daq_objs, channel=channel, device_id=device_id)
            time.sleep(wait_time)

    return daq_objs

def initialize_zurich_lockin(device_id=DEVICE_ID):
    apilevel=6
    objs = ziutils.create_api_session(device_id, apilevel)
    daq, device, props = objs

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


##################
### Scope read ###
##################

def read_zurich_spectrum(daq_objs=None):
    '''
    reads scope on lockin and return dictionary with scope data. In this primitive form, up to user to make sure the that the scope has been setup correctly.

    UNDER CONSTRUCTION
    '''

    # initialize
    if daq_objs is None:
        daq, device, props = initialize_zurich_lockin()
        daq_objs = [daq, device, props]
    else:
        daq, device, props = daq_objs

    # subscribe to channels and read mfli
    for channel in channels:
        daq.subscribe(f'/{device}/scopes/{channel-1}/sample')

    # read lockin
    daq.sync()
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
        autophase = daq.getDouble(f'/{device}/demods/{channel-1}/phaseshift')
        osc = daq.getInt(f'/{device}/demods/{channel-1}/oscselect')
        freq = daq.getDouble(f'/{device}/oscs/{osc}/freq')
        data_dict[f'Demod {channel} x'] = x
        data_dict[f'Demod {channel} y'] = y
        data_dict[f'Demod {channel} r'] = np.abs(x + 1j*y)
        data_dict[f'Demod {channel} phase'] = np.arctan2(y,x)
        data_dict[f'Demod {channel} autophase'] = autophase
        data_dict[f'Demod {channel} oscillator'] = osc+1
        data_dict[f'Demod {channel} frequency'] = freq

    # add any other data from lockin to data_dict

    return data_dict

#########################
### Wrapper functions ###
#########################

###
### Lockin 1
###
def read_zurich_lockin1(daq_objs=None, time_constant=0.3, poll_length=None, poll_timeout=500, wait_factor=3, channels=[1,2,3,4], name='LI1'):
    return read_zurich_lockin(daq_objs=daq_objs, time_constant=time_constant, poll_length=poll_length, poll_timeout=poll_timeout, wait_factor=wait_factor, channels=channels, device_id=DEVICE_ID, name=name)

def set_zurich_aux_offset_1(offset, daq_objs=None, wait_time=0, check_stability=True):
    return set_zurich_aux_offset(offset, daq_objs, wait_time, 1)

def set_zurich_aux_offset_2(offset, daq_objs=None, wait_time=0):
    set_zurich_aux_offset(offset, daq_objs, wait_time, 2)

def get_zurich_aux_offset_1(daq_objs=None, wait_time=0):
    return get_zurich_aux_offset(daq_objs, wait_time, 1)

def get_zurich_aux_offset_2(daq_objs=None, wait_time=0):
    return get_zurich_aux_offset(daq_objs, wait_time, 2)

def initialize_zurich_lockin1():
    return initialize_zurich_lockin(DEVICE_ID)

###
### Lockin 2
###

def read_zurich_lockin2(daq_objs=None, time_constant=0.3, poll_length=None, poll_timeout=500, wait_factor=3, channels=[1], name='LI2'):
    return read_zurich_lockin(daq_objs=daq_objs, time_constant=time_constant, poll_length=poll_length, poll_timeout=poll_timeout, wait_factor=wait_factor, channels=channels, device_id=DEVICE_ID2, name=name)

def set_zurich_output_amplitude2(amplitude, daq_objs=None, wait_time=0, channel=2, check_stability=True):
    return set_zurich_output_amplitude(amplitude, daq_objs=daq_objs, wait_time=wait_time, channel=channel, check_stability=check_stability, device_id=DEVICE_ID2)

def set_zurich_output_ramp2(v, daq_objs=None, vstep=0.5, channel=2, wait_time=2, check_stability=True):
    return set_zurich_output_ramp(v, daq_objs=daq_objs, vstep=vstep, channel=channel, wait_time=wait_time, check_stability=check_stability, device_id=DEVICE_ID2)

def read_zurich_output_amplitude2(daq_objs=None, channel=2):
    return read_zurich_output_amplitude(daq_objs=daq_objs, channel=channel, device_id=DEVICE_ID2)

def set_zurich_osc2(osc, daq_objs=None, wait_time=0, channel=1, check_stability=True):
    return set_zurich_osc(osc, daq_objs=daq_objs, wait_time=wait_time, channel=channel, check_stability=True, device_id=DEVICE_ID2)

def read_zurich_osc2(daq_objs=None, channel=1):
    return read_zurich_osc(daq_objs=daq_objs, channel=channel, device_id=DEVICE_ID2)

def set_zurich_frequency2(freq, daq_objs=None, wait_time=0, osc=1):
    return set_zurich_frequency(freq, daq_objs=daq_objs, wait_time=wait_time, osc=osc, device_id=DEVICE_ID2)

def read_zurich_frequency2(daq_objs=None, osc=1):
    return read_zurich_frequency(daq_objs=daq_objs, osc=osc, device_id=DEVICE_ID2)

def initialize_zurich_lockin2():
    return initialize_zurich_lockin(DEVICE_ID2)
