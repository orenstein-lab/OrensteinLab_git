from OrensteinLab_git.Measurement.Alex.control import attocube, lakeshore, newport, opticool, razorbill, zurich, montana_temp

'''
This file contains system specific motor and instrument configurations.

    - motor_dict: a dictionary for defining system motors, which form the backbone of an experiment

    - insrtument_dict: a dictionary for defining system measurement instruments, such as a lockin or spectrometer

    - default_record_dict: a dictionary for motors that should be automatically recorded and added to each file
'''

#####################
### System Motors ###
#####################
# entries of the form motor:{'move':move_function, 'read':read_function, 'init':initialize_function, 'close':close_function}
motor_dict = {
'x':{'move':attocube.move_x, 'read':attocube.read_x, 'init':attocube.initialize_attocube, 'close':attocube.close_attocube, 'move_back':10, 'name':'x (um)'},

'y':{'move':attocube.move_y, 'read':attocube.read_y, 'init':attocube.initialize_attocube, 'close':attocube.close_attocube, 'move_back':10, 'name':'y (um)'},

'z':{'move':attocube.move_z, 'read':attocube.read_z, 'init':attocube.initialize_attocube, 'close':attocube.close_attocube, 'move_back':10, 'name':'z (um)'},

'temp':{'move':lakeshore.set_temperature, 'read':lakeshore.read_temperature, 'init':lakeshore.initialize_lakeshore, 'close':lakeshore.close_lakeshore, 'move_back':0, 'name':'Temperature (K)'},

'axis_1':{'move':newport.move_axis_1, 'read':newport.read_axis_1, 'init':newport.initialize_axis_1, 'close':newport.close_axis_1, 'move_back':1, 'name':'Angle 1 (deg)'},

'axis_2':{'move':newport.move_axis_2, 'read':newport.read_axis_2, 'init':newport.initialize_axis_2, 'close':newport.close_axis_2, 'move_back':1, 'name':'Angle 2 (deg)'},

'axis_3':{'move':newport.move_axis_3, 'read':newport.read_axis_3, 'init':newport.initialize_axis_3, 'close':newport.close_axis_3, 'move_back':1, 'name':'Angle 3 (deg)'},

'delay_stage':{'move':newport.move_axis_3, 'read':newport.read_axis_3, 'init':newport.initialize_axis_3, 'close':newport.close_axis_3, 'move_back':1, 'name':'Delay Stage (mm)'},

'strain_cap':{'move':razorbill.set_strain_capacitance, 'read':razorbill.read_strain_capacitance, 'init':razorbill.initialize_strain_cell_client, 'close':razorbill.close_strain_cell_client, 'move_back':0, 'name':'Capacitance (pF)'},

'strain_ps':{'move':razorbill.set_strain_ps, 'read':razorbill.read_strain_ps, 'init':razorbill.initialize_strain_cell_client, 'close':razorbill.close_strain_cell_client, 'move_back':0, 'name':'Voltage (V)'},

'zurich_output':{'move':zurich.set_zurich_output_amplitude, 'read':zurich.read_zurich_output_amplitude, 'init':zurich.initialize_zurich_lockin, 'close':zurich.close_zurich_lockin, 'move_back':0, 'name':'Lock-in Ouput Voltage (V)'},

'zurich_frequency':{'move':zurich.set_zurich_frequency, 'read':zurich.read_zurich_frequency, 'init':zurich.initialize_zurich_lockin, 'close':zurich.close_zurich_lockin, 'move_back':0, 'name':'Lock-in Frequency (Hz)'},

'opticool_temp':{'move':opticool.set_opticool_temperature, 'read':opticool.read_opticool_temperature, 'init':opticool.initialize_opticool, 'close':opticool.close_opticool, 'move_back':0, 'name':'Opticool Temperature (K)'},

'opticool_field':{'move':opticool.set_opticool_field, 'read':opticool.read_opticool_field, 'init':opticool.initialize_opticool, 'close':opticool.close_opticool, 'move_back':0, 'name':'Field (Oe)'},

'corotate_axes12':{'move':newport.corotate_axes12, 'read':newport.read_corotate_axes12, 'init':newport.initialize_corotate_axes12, 'close':newport.close_corotate_axes12, 'move_back':1, 'name':'Corotation Axes (deg)'},

'corotate_axes_unbal':{'move':newport.corotate_axes12, 'read':newport.read_corotate_axes12, 'init':newport.initialize_corotate_axes12, 'close':newport.close_corotate_axes12, 'move_back':1, 'name':'Unbalanced Corotation (deg)'},

'galvo_x':{'move':zurich.set_zurich_aux_offset_1, 'read':zurich.get_zurich_aux_offset_1, 'init':zurich.initialize_zurich_lockin, 'close':zurich.close_zurich_lockin, 'move_back':0, 'name':'Galvo x Voltage (V)'},

'galvo_y':{'move':zurich.set_zurich_aux_offset_2, 'read':zurich.get_zurich_aux_offset_2, 'init':zurich.initialize_zurich_lockin, 'close':zurich.close_zurich_lockin, 'move_back':0, 'name':'Galvo y Voltage (V)'},

'montana_temp':{'move':montana_temp.set_montana_temperature, 'read':montana_temp.read_montana_temperature, 'init':montana_temp.initialize_montana, 'close':montana_temp.close_montana, 'move_back':0, 'name':'Montana Platform Temperature (K)'}
}

meta_motors = ['x', 'y', 'z', 'axis_1', 'axis_2', 'axis_3', 'temp', 'zurich_frequency', 'strain_ps', 'strain_cap', 'montana_temp']
#meta_motors=[]

##########################
### System Instruments ###
##########################
# entries of the form motor:{'read':read_function, 'init':initialize_functions}
instrument_dict = {
'zurich_lockin':{'read':zurich.read_zurich_lockin, 'init':zurich.initialize_zurich_lockin}
}

default_record_dict = {}
