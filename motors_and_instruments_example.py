'''
This file contains system specific motor and instrument configurations.

    - MOTOR_DICT: a dictionary for defining system motors, which form the backbone of an experiment

    - insrtument_dict: a dictionary for defining system measurement instruments, such as a lockin or spectrometer

    - META_MOTORS: motors to record into metadata header

    - ACTIVE_MOTORS: motors to measure when acquiring data

    - DEFAULT_VARS: names of variable to plot during scans by default
'''

############################
### Import Devices Files ###
############################

#####################
### System Motors ###
#####################
# entries of the form motor:{'move':move_function, 'read':read_function, 'init':initialize_function, 'close':close_function, 'move_back':0, 'name':motor_name}
MOTOR_DICT = {}

##########################
### System Instruments ###
##########################
# entries of the form instrument:{'read':read_function, 'init':initialize_functions, 'close':close_function, 'name':instrument_name}
INSTRUMENT_DICT = {}

# keys in MOTOR_DICT for motors to measure during measurements
ACTIVE_MOTORS = []

# instruments to use during measurement
ACTIVE_INSTRUMENTS = []

# names of variables to display during plotting by default
DEFAULT_VARS = []
