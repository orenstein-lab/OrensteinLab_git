'''
Wrapper functions which recreate old versions of measurement functions (ie, with same arguments, behavior of functions may differ slightly)
'''
get_ipython().run_line_magic('matplotlib', 'notebook')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
import time
import threading
from tqdm.auto import tqdm
import scipy.optimize as opt
import scipy.interpolate as interp
import inspect
import pickle
import OrensteinLab_git.devices.newport as newport
from OrensteinLab_git.devices.concurrency_classes import LockedVar, StoppableThread, LockedDict
from OrensteinLab_git.motors_and_instruments import MOTOR_DICT, INSTRUMENT_DICT, META_MOTORS, ACTIVE_MOTORS, ACTIVE_INSTRUMENTS, DEFAULT_VARS
import OrensteinLab_git.helper as helper
import OrensteinLab_git.plotters as plotters
from orenstein_analysis.measurement import loader, process
from OrensteinLab_git import measurement as meas

################
### Medthods ###
################

def lockin_time_series(recording_time, filename_head=None, filename=None, mobj_measure_dict={}, metadata={}, time_constant=0.3, channel_index=1, savefile=True, daq_objs=None, plot_flag=True, override_metadata=False):
    inst_dict={'zurich_lockin':{'time_constant':time_constant, 'channel_index':1}}
    iobj_dict={'zurich_lockin':daq_objs}
    vars = [f'Demod {channel_index} x', f'Demod {channel_index} y', f'Demod {channel_index} r']
    close_devices=False
    if mobj_measure_dict=={}:
        close_devices=True
    meas.timed_measurement(recording_time, inst_dict=inst_dict, mobj_dict=mobj_measure_dict, iobj_dict=iobj_dict, filename_head=filename_head, filename=filename, savefile=savefile, plot=plot_flag, vars=vars, metadata=metadata, close_devices=close_devices)
    
def rotate_scan(start_angle, end_angle, step_size, filename_head=None, filename=None, axis_index=1, mobj_measure_dict={}, metadata={}, showplot=True, time_constant=0.3, channel_index=1, R_channel_index=4, daq_objs=None, axis_1=None, axis_2=None, savefile=True, override_metadata=False):
    inst_dict={'zurich_lockin':{'time_constant':time_constant, 'channel_index':1}}
    iobj_dict={'zurich_lockin':daq_objs}
    vars = [f'Demod {channel_index} x', f'Demod {channel_index} y', f'Demod {channel_index} r']
    close_devices=False
    if mobj_measure_dict=={}:
        close_devices=True
    meas.rotate_scan(start_angle, end_angle, step_size, axis_index=axis_index, inst_dict=inst_dict, mobj_dict=mobj_measure_dict, iobj_dict=iobj_dict, filename_head=filename_head, filename=filename, savefile=savefile, print_flag=False, vars=vars, plot=showplot, metadata=metadata, close_devices=close_devices)

def corotate_scan(start_angle, end_angle, step_size, angle_offset, rate_axis_2=1, filename_head=None, filename=None, mobj_measure_dict={}, metadata={}, showplot=True, time_constant=0.3, channel_index=1, R_channel_index=4, daq_objs=None, axis_1=None, axis_2=None, savefile=True, override_metadata=False):
    inst_dict={'zurich_lockin':{'time_constant':time_constant, 'channel_index':1}}
    iobj_dict={'zurich_lockin':daq_objs}
    vars = [f'Demod {channel_index} x', f'Demod {channel_index} y', f'Demod {channel_index} r']
    close_devices=False
    if mobj_measure_dict=={}:
        close_devices=True
    meas.corotate_scan(start_angle, end_angle, step_size, angle_offset, rate_axis_2=rate_axis_2, inst_dict=inst_dict, mobj_dict=mobj_measure_dict, iobj_dict=iobj_dict, filename_head=filename_head, filename=filename, savefile=savefile, print_flag=False, vars=vars, plot=showplot, metadata=metadata, close_devices=close_devices)

def motor_scan(map_dict, filename_head=None, filename=None, showplot=True, time_constant=0.3, channel_index=1, R_channel_index=2, print_flag=False, savefile=True, metadata={}, daq_objs=None):
    inst_dict={'zurich_lockin':{'time_constant':time_constant, 'channel_index':1}}
    iobj_dict={'zurich_lockin':daq_objs}
    vars = [f'Demod {channel_index} x', f'Demod {channel_index} y', f'Demod {channel_index} r']
    meas.motor_scan(map_dict, inst_dict=inst_dict, mobj_dict={}, iobj_dict=iobj_dict, filename_head=filename_head, filename=filename, savefile=savefile, print_flag=print_flag, vars=vars, plot=showplot, metadata=metadata, close_devices=True)

def motor_scan_balance(map_dict, balance, balance_table=None, slope=0, tol=0, balance_channel=3, autobalance_flag=True, filename_head=None, filename=None, showplot=True, time_constant=0.3, channel_index=1, R_channel_index=4, print_flag=False, savefile=True, metadata={}, daq_objs=None):
    # under construction
    return 0

def rotate_map(map_dict, start_angle, end_angle, step_size, filename_head=None, filename=None, axis_index=1, showplot=False, time_constant=0.3, channel_index=1, R_channel_index=4, daq_objs=None, print_flag=False, savefile=True, metadata={}):
    inst_dict={'zurich_lockin':{'time_constant':time_constant, 'channel_index':1}}
    iobj_dict={'zurich_lockin':daq_objs}
    vars = [f'Demod {channel_index} x', f'Demod {channel_index} y', f'Demod {channel_index} r']
    meas.rotate_map(map_dict, start_angle, end_angle, step_size, channel_index=channel_index, inst_dict=inst_dict, mobj_dict={}, iobj_dict=iobj_dict, filename_head=filename_head, filename=filename, savefile=savefile, print_flag=print_flag, vars=vars, plot=showplot, metadata=metadata, close_devices=True)

def corotate_map(map_dict, start_angle, end_angle, step_size, angle_offset, rate_axis_2=1, filename_head=None, filename=None, measure_motors=[], showplot=False, time_constant=0.3, channel_index=1, R_channel_index=4, print_flag=False, savefile=True, metadata={}):
    inst_dict={'zurich_lockin':{'time_constant':time_constant, 'channel_index':1}}
    iobj_dict={'zurich_lockin':daq_objs}
    vars = [f'Demod {channel_index} x', f'Demod {channel_index} y', f'Demod {channel_index} r']
    meas.corotate_map(map_dict, start_angle, end_angle, step_size, angle_offset, rate_axis_2=1, inst_dict={}, mobj_dict={}, iobj_dict={}, filename_head=None, filename=None, savefile=True, print_flag=False, vars=[], plot=True, metadata={}, close_devices=True)

def cont_motor_meas(sweep_list,filename_head=None, filename=None, mobj_measure_dict={}, metadata={}, time_constant=0.3, channel_index=1, savefile=True, daq_objs=None, plot_flag=True, override_metadata=False):
    inst_dict={'zurich_lockin':{'time_constant':time_constant, 'channel_index':1}}
    iobj_dict={'zurich_lockin':daq_objs}
    vars = [f'Demod {channel_index} x', f'Demod {channel_index} y', f'Demod {channel_index} r']
    close_devices=False
    if mobj_measure_dict=={}:
        close_devices=True
    meas.motor_sequence(sequence_list, inst_dict=inst_dict, mobj_dict=mobj_measure_dict, iobj_dict=iobj_dict, filename_head=filename_head, filename=filename, savefile=savefile, plot=showplot, vars=vars, metadata=metadata, close_devices=close_devices)

#########################
### Balancing Methods ###
#########################

def find_balance_angle(start_angle, end_angle, step_size, balance_at=0, offset=0, go_to_balance_angle=True, axis_index=2, channel_index=1, time_constant=0.3, R_channel_index=2, daq_objs=None, axis_1=None, axis_2=None):
    inst_dict={'zurich_lockin':{'time_constant':time_constant, 'channel_index':1}}
    iobj_dict={'zurich_lockin':daq_objs}
    mobj_dict = {'axis_1':axis_1, 'axis_2':axis_2}
    var = f'Demod {channel_index} x'
    meas.find_balance_angle(start_angle, end_angle, step_size, balance_at=balance_at, offset=offset, bal_axis=axis_index, go_to_balance_angle=go_to_balance_angle, var=var, inst_dict=inst_dict, mobj_dict=mobj_dict, iobj_dict=iobj_dict)

def autobalance(slope, tolerance, daq_objs=None, axis_1=None, axis_2=None, balance_at=None, offset=0, channel_index=1, time_constant=0.3, print_flag=True):
    inst_dict={'zurich_lockin':{'time_constant':time_constant, 'channel_index':1}}
    iobj_dict={'zurich_lockin':daq_objs}
    mobj_dict = {'axis_1':axis_1, 'axis_2':axis_2}
    var = f'Demod {channel_index} x'
    meas.autobalance(slope, tolerance, offset=offset, bal_axis=2, var=var, inst_dict=inst_dict, mobj_dict=mobj_dict, iobj_dict=iobj_dict, print_flag=print_flag)