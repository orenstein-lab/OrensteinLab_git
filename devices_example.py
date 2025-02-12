'''
This file contains system specific motor and instrument configurations.

    - motor_dict: a dictionary for defining system motors, which form the backbone of an experiment

    - insrtument_dict: a dictionary for defining system measurement instruments, such as a lockin or spectrometer

    - meta_motors: motors to record into metadata header

    - active_motors: motors to measure when acquiring data
'''

#####################
### System Motors ###
#####################
# entries of the form motor:{'move':move_function, 'read':read_function, 'init':initialize_function, 'close':close_function, 'move_back':0, 'name':motor_name}
motor_dict = {}

##########################
### System Instruments ###
##########################
# entries of the form instrument:{'read':read_function, 'init':initialize_functions, 'close':close_function, 'name':instrument_name}
instrument_dict = {}

# keys in motor_dict for motors to be recorded into metadata header - useful if there are motors which you don't want to constantly stream but do want a record of
meta_motors = []

# keys in motor_dict for motors to measure during measurements
active_motors = []

# instruments to use during measurement
active_instruments = []

# names of variables to return by default during scans and to display during plotting
default_vars = []

