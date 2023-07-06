from strain_control.strain_client import StrainClient
import OrensteinLab_git.Measurement.Alex.control as ctrl
import OrensteinLab_git.Instrument.montana.cryocore as cryocore
import os

'''
This file contains system specific configurations.
s
    - config_dict: a dictionary containing any configurations the user may want to specifiy, such as hardware ports for specific instruments

    - motor_dict: a dictionary for defining system motors, which form the backbone of an experiment

    - insrtument_dict: a dictionary for defining system measurement instruments, such as a lockin or spectrometer

    - default_record_dict: a dictionary for motors that should be automatically recorded and added to each file
'''

############################
### System Configuration ###
############################
config_dict = {
'Zurich Lockin ID':'dev3425',
'ESP COM PORT':'COM4',
'ESP Model':301,
'Lakeshore Model':336,
'Attocube Handle': f'{os.path.dirname(__file__)}\\attocube_handle',
'Opticool IP':'131.243.163.240',
'Opticool Port':'5000'
}

#####################
### System Motors ###
#####################
# entries of the form motor:{'move':move_function, 'read':read_function, 'init':initialize_function, 'close':close_function}
motor_dict = {
'x':{'move':ctrl.move_x, 'read':ctrl.read_x, 'init':ctrl.initialize_attocube, 'close':ctrl.close_attocube, 'move_back':10, 'name':'x (um)'},

'y':{'move':ctrl.move_y, 'read':ctrl.read_y, 'init':ctrl.initialize_attocube, 'close':ctrl.close_attocube, 'move_back':10, 'name':'y (um)'},

'z':{'move':ctrl.move_z, 'read':ctrl.read_z, 'init':ctrl.initialize_attocube, 'close':ctrl.close_attocube, 'move_back':10, 'name':'z (um)'},

'temp':{'move':ctrl.set_temperature, 'read':ctrl.read_temperature, 'init':ctrl.initialize_lakeshore, 'close':ctrl.close_lakeshore, 'move_back':0, 'name':'Temperature (K)'},

'axis_1':{'move':ctrl.rotate_axis_1, 'read':ctrl.read_axis_1, 'init':ctrl.initialize_rot_axis_1, 'close':ctrl.close_rot_axis_1, 'move_back':1, 'name':'Angle 1 (deg)'},

'axis_2':{'move':ctrl.rotate_axis_2, 'read':ctrl.read_axis_2, 'init':ctrl.initialize_rot_axis_2, 'close':ctrl.close_rot_axis_2, 'move_back':1, 'name':'Angle 2 (deg)'},

'axis_3':{'move':ctrl.rotate_axis_3, 'read':ctrl.read_axis_3, 'init':ctrl.initialize_rot_axis_3, 'close':ctrl.close_rot_axis_3, 'move_back':1, 'name':'Angle 2 (deg)'},

'delay_stage':{'move':ctrl.rotate_axis_3, 'read':ctrl.read_axis_3, 'init':ctrl.initialize_rot_axis_3, 'close':ctrl.close_rot_axis_3, 'move_back':1, 'name':'Delay Stage (mm)'},

'strain_cap':{'move':ctrl.set_strain_capacitance, 'read':ctrl.read_strain_capacitance, 'init':ctrl.initialize_strain_cell_client, 'close':ctrl.close_strain_cell_client, 'move_back':0, 'name':'Capacitance (pF)'},

'strain_ps':{'move':ctrl.set_strain_ps, 'read':ctrl.read_strain_ps, 'init':ctrl.initialize_strain_cell_client, 'close':ctrl.close_strain_cell_client, 'move_back':0, 'name':'Voltage (V)'},

'zurich_output':{'move':ctrl.set_zurich_output_amplitude, 'read':ctrl.read_zurich_output_amplitude, 'init':ctrl.initialize_zurich_lockin, 'close':ctrl.close_zurich_lockin, 'move_back':0, 'name':'Lock-in Ouput Voltage (V)'},

'zurich_frequency':{'move':ctrl.set_zurich_frequency, 'read':ctrl.read_zurich_frequency, 'init':ctrl.initialize_zurich_lockin, 'close':ctrl.close_zurich_lockin, 'move_back':0, 'name':'Lock-in Frequency (Hz)'},

'opticool_temp':{'move':ctrl.set_opticool_temperature, 'read':ctrl.read_opticool_temperature, 'init':ctrl.initialize_opticool, 'close':ctrl.close_opticool, 'move_back':0, 'name':'Temperature (K)'},

'opticool_field':{'move':ctrl.set_opticool_field, 'read':ctrl.read_opticool_field, 'init':ctrl.initialize_opticool, 'close':ctrl.close_opticool, 'move_back':0, 'name':'Field (Oe)'},

'corotate_axes12':{'move':ctrl.corotate_axes12, 'read':ctrl.read_corotate_axes12, 'init':ctrl.initialize_corotate_axes12, 'close':ctrl.close_corotate_axes12, 'move_back':1, 'name':'Corotation Axes (deg)'}
}

##########################
### System Instruments ###
##########################
# entries of the form motor:{'read':read_function, 'init':initialize_functions}
instrument_dict = {
'zurich_lockin':{'read':ctrl.read_zurich_lockin, 'init':ctrl.initialize_zurich_lockin}
}

default_record_dict = {}
