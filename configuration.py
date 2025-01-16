import os

'''
This file contains system specific configurations.
s
    - config_dict: a dictionary containing any configurations the user may want to specifiy, such as hardware ports for specific instruments

'''

############################
### System Configuration ###
############################
config_dict = {
'Zurich Lockin ID':'dev3232',
'ESP COM Port':'COM4',
'ESP Model':'301',
'Lakeshore Model':'335',
'Attocube Handle': f'{os.path.dirname(__file__)}\\attocube_handle',
'AMC300 IP': "10.1.1.22",
'Zurich Handle': f'{os.path.dirname(__file__)}\\zurich_handle',
'Opticool Handle': f'{os.path.dirname(__file__)}\\opticool_handle',
'Opticool IP':'131.243.163.240',
'Opticool Port':'5000',
'Montana Address':'10.1.1.15'
}
